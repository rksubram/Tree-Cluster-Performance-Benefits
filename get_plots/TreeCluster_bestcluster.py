#!/usr/bin/env python3
from scipy.optimize import minimize
import sys
import numpy as np
from math import log
from niemads import DisjointSet
from queue import PriorityQueue,Queue
from treeswift import read_tree_newick
from sys import stderr
NUM_THRESH = 1000 # number of thresholds for the threshold-free methods to use

# merge two sorted lists into a sorted list
def merge_two_sorted_lists(x,y):
    out = list(); i = 0; j = 0
    while i < len(x) and j < len(y):
        if x[i] < y[j]:
            out.append(x[i]); i+= 1
        else:
            out.append(y[j]); j += 1
    while i < len(x):
        out.append(x[i]); i += 1
    while j < len(y):
        out.append(y[j]); j += 1
    return out

# merge multiple sorted lists into a sorted list
def merge_multi_sorted_lists(lists):
    pq = PriorityQueue()
    for l in range(len(lists)):
        if len(lists[l]) != 0:
            pq.put((lists[l][0],l))
    inds = [1 for _ in range(len(lists))]
    out = list()
    while not pq.empty():
        d,l = pq.get(); out.append(d)
        if inds[l] < len(lists[l]):
            pq.put((lists[l][inds[l]],l)); l += 1
    return out

# get the median of a sorted list
def median(x):
    if len(x) % 2 != 0:
        return x[int(len(x)/2)]
    else:
        return (x[int(len(x)/2)]+x[int(len(x)/2)-1])/2

# get the average of a list
def avg(x):
    return float(sum(x))/len(x)

# convert p-distance to Jukes-Cantor distance
def p_to_jc(d,seq_type):
    b = {'dna':3./4., 'protein':19./20.}[seq_type]
    return -1*b*log(1-(d/b))

# cut out the current node's subtree (by setting all nodes' DELETED to True) and return list of leaves
def cut(node):
    cluster = list()
    descendants = Queue(); descendants.put(node)
    while not descendants.empty():
        descendant = descendants.get()
        if descendant.DELETED:
            continue
        descendant.DELETED = True
        descendant.left_dist = 0; descendant.right_dist = 0; descendant.edge_length = 0
        if descendant.is_leaf():
            cluster.append(str(descendant))
        else:
            for c in descendant.children:
                descendants.put(c)
    return cluster

# initialize properties of input tree and return set containing taxa of leaves
def prep(tree,support):
    tree.resolve_polytomies(); tree.suppress_unifurcations()
    leaves = set()
    for node in tree.traverse_postorder():
        if node.edge_length is None:
            node.edge_length = 0
        node.DELETED = False
        if node.is_leaf():
            leaves.add(str(node))
        else:
            try:
                node.confidence = float(str(node))
            except:
                node.confidence = 100. # give edges without support values support 100
            if node.confidence < support: # don't allow low-support edges
                node.edge_length = float('inf')
    return leaves

# return a sorted list of all unique pairwise leaf distances <= a given threshold
def pairwise_dists_below_thresh(tree,threshold):
    pairwise_dists = set()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.leaf_dists = {0}; node.min_leaf_dist = 0
        else:
            children = list(node.children)
            for i in range(len(children)-1):
                c1 = children[i]
                for j in range(i+1,len(children)):
                    c2 = children[j]
                    for d1 in c1.leaf_dists:
                        for d2 in c2.leaf_dists:
                            pd = d1 + c1.edge_length + d2 + c2.edge_length
                            if pd <= threshold:
                                pairwise_dists.add(pd)
            node.leaf_dists = set(); node.min_leaf_dist = float('inf')
            for c in children:
                if c.min_leaf_dist + c.edge_length > threshold:
                    continue
                for d in c.leaf_dists:
                    nd = d+c.edge_length
                    if nd < threshold:
                        node.leaf_dists.add(nd)
                    if nd < node.min_leaf_dist:
                        node.min_leaf_dist = nd
    return sorted(pairwise_dists)

# split leaves into minimum number of clusters such that the maximum leaf pairwise distance is below some threshold
def min_clusters_threshold_max(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()
    for node in tree.traverse_postorder():
        # if I've already been handled, ignore me
        if node.DELETED:
            continue

        # find my undeleted max distances to leaf
        if node.is_leaf():
            node.left_dist = 0; node.right_dist = 0
        else:
            children = list(node.children)
            if children[0].DELETED and children[1].DELETED:
                cut(node); continue
            if children[0].DELETED:
                node.left_dist = 0
            else:
                node.left_dist = max(children[0].left_dist,children[0].right_dist) + children[0].edge_length
            if children[1].DELETED:
                node.right_dist = 0
            else:
                node.right_dist = max(children[1].left_dist,children[1].right_dist) + children[1].edge_length

            # if my kids are screwing things up, cut out the longer one
            if node.left_dist + node.right_dist > threshold:
                if node.left_dist > node.right_dist:
                    cluster = cut(children[0])
                    node.left_dist = 0
                else:
                    cluster = cut(children[1])
                    node.right_dist = 0

                # add cluster
                if len(cluster) != 0:
                    clusters.append(cluster)
                    for leaf in cluster:
                        leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# median leaf pairwise distance cannot exceed threshold, and clusters must define clades
def min_clusters_threshold_med_clade(tree,threshold,support):
    leaves = prep(tree,support)
    # bottom-up traversal to compute median pairwise distances
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.med_pair_dist = 0
            node.leaf_dists = [0]
            node.pair_dists = list()
        else:
            children = list(node.children)
            l_leaf_dists = [d + children[0].edge_length for d in children[0].leaf_dists]
            r_leaf_dists = [d + children[1].edge_length for d in children[1].leaf_dists]
            node.leaf_dists = merge_two_sorted_lists(l_leaf_dists,r_leaf_dists)
            if len(l_leaf_dists) < len(r_leaf_dists):
                across_leaf_dists = [[l+r for r in r_leaf_dists] for l in l_leaf_dists]
            else:
                across_leaf_dists = [[l+r for l in l_leaf_dists] for r in r_leaf_dists]
            node.pair_dists = merge_multi_sorted_lists([children[0].pair_dists,children[1].pair_dists] + across_leaf_dists)
            if node.pair_dists[-1] == float('inf'):
                node.med_pair_dist = float('inf')
            else:
                node.med_pair_dist = median(node.pair_dists)
            for c in (children[0],children[1]):
                del c.leaf_dists; del c.pair_dists

    # top-down traversal to cut out clusters
    clusters = list()
    traverse = Queue(); traverse.put(tree.root)
    while not traverse.empty():
        node = traverse.get()
        if node.med_pair_dist <= threshold:
            clusters.append(cut(node))
        else:
            for c in node.children:
                traverse.put(c)
    return clusters

# average leaf pairwise distance cannot exceed threshold, and clusters must define clades
def min_clusters_threshold_avg_clade(tree,threshold,support):
    leaves = prep(tree,support)
    # bottom-up traversal to compute average pairwise distances
    for node in tree.traverse_postorder():
        node.total_pair_dist = 0; node.total_leaf_dist = 0
        if node.is_leaf():
            node.num_leaves = 1
            node.avg_pair_dist = 0
        else:
            children = list(node.children)
            node.num_leaves = sum(c.num_leaves for c in children)
            node.total_pair_dist = children[0].total_pair_dist + children[1].total_pair_dist + (children[0].total_leaf_dist*children[1].num_leaves + children[1].total_leaf_dist*children[0].num_leaves)
            node.total_leaf_dist = (children[0].total_leaf_dist + children[0].edge_length*children[0].num_leaves) + (children[1].total_leaf_dist + children[1].edge_length*children[1].num_leaves)
            node.avg_pair_dist = node.total_pair_dist/((node.num_leaves*(node.num_leaves-1))/2)

    # top-down traversal to cut out clusters
    clusters = list()
    traverse = Queue(); traverse.put(tree.root)
    while not traverse.empty():
        node = traverse.get()
        if node.avg_pair_dist <= threshold:
            clusters.append(cut(node))
        else:
            for c in node.children:
                traverse.put(c)
    return clusters

# total branch length cannot exceed threshold, and clusters must define clades
def min_clusters_threshold_sum_bl_clade(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.left_total = 0; node.right_total = 0
        else:
            children = list(node.children)
            if children[0].DELETED and children[1].DELETED:
                cut(node); continue
            if children[0].DELETED:
                node.left_total = 0
            else:
                node.left_total = children[0].left_total + children[0].right_total + children[0].edge_length
            if children[1].DELETED:
                node.right_total = 0
            else:
                node.right_total = children[1].left_total + children[1].right_total + children[1].edge_length
            if node.left_total + node.right_total > threshold:
                cluster_l = cut(children[0])
                node.left_total = 0
                cluster_r = cut(children[1])
                node.right_total = 0
                for cluster in (cluster_l,cluster_r):
                    if len(cluster) != 0:
                        clusters.append(cluster)
                        for leaf in cluster:
                            leaves.remove(leaf)
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# total branch length cannot exceed threshold
def min_clusters_threshold_sum_bl(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.left_total = 0; node.right_total = 0
        else:
            children = list(node.children)
            if children[0].DELETED and children[1].DELETED:
                cut(node); continue
            if children[0].DELETED:
                node.left_total = 0
            else:
                node.left_total = children[0].left_total + children[0].right_total + children[0].edge_length
            if children[1].DELETED:
                node.right_total = 0
            else:
                node.right_total = children[1].left_total + children[1].right_total + children[1].edge_length
            if node.left_total + node.right_total > threshold:
                if node.left_total > node.right_total:
                    cluster = cut(children[0])
                    node.left_total = 0
                else:
                    cluster = cut(children[1])
                    node.right_total = 0
                if len(cluster) != 0:
                    clusters.append(cluster)
                    for leaf in cluster:
                        leaves.remove(leaf)
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# single-linkage clustering using Metin's cut algorithm
def single_linkage_cut(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()

	# find closest leaf below (dist,leaf)
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.min_below = (0,node.label)
        else:
            node.min_below = min((c.min_below[0]+c.edge_length,c.min_below[1]) for c in node.children)

    # find closest leaf above (dist,leaf)
    for node in tree.traverse_preorder():
        node.min_above = (float('inf'),None)
        if node.is_root():
            continue
        # min distance through sibling
        for c in node.parent.children:
            if c != node:
                dist = node.edge_length + c.edge_length + c.min_below[0]
                if dist < node.min_above[0]:
                    node.min_above = (dist,c.min_below[1])
        # min distance through grandparent
        if not c.parent.is_root():
            dist = node.edge_length + node.parent.min_above[0]
            if dist < node.min_above[0]:
                node.min_above = (dist,node.parent.min_above[1])

    # find clusters
    for node in tree.traverse_postorder(leaves=False):
        # assume binary tree here (prep function guarantees this)
        l_child,r_child = node.children
        l_dist = l_child.min_below[0] + l_child.edge_length
        r_dist = r_child.min_below[0] + r_child.edge_length
        a_dist = node.min_above[0]
        bad = [0,0,0] # left, right, up
        if l_dist + r_dist > threshold:
            bad[0] += 1; bad[1] += 1
        if l_dist + a_dist > threshold:
            bad[0] += 1; bad[2] += 1
        if r_dist + a_dist > threshold:
            bad[1] += 1; bad[2] += 1
        # cut either (or both) children
        for i in [0,1]:
            if bad[i] == 2:
                cluster = cut(node.children[i])
                if len(cluster) != 0:
                    clusters.append(cluster)
                    for leaf in cluster:
                        leaves.remove(leaf)
        # cut above (equals cutting me)
        if bad[2] == 2: # if cutting above, just cut me
            cluster = cut(node)
            if len(cluster) != 0:
                clusters.append(cluster)
                for leaf in cluster:
                    leaves.remove(leaf)
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# single-linkage clustering using Niema's union algorithm
def single_linkage_union(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()

    # find closest leaf below (dist,leaf)
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.min_below = (0,node.label)
        else:
            node.min_below = min((c.min_below[0]+c.edge_length,c.min_below[1]) for c in node.children)

    # find closest leaf above (dist,leaf)
    for node in tree.traverse_preorder():
        node.min_above = (float('inf'),None)
        if node.is_root():
            continue
        # min distance through sibling
        for c in node.parent.children:
            if c != node:
                dist = node.edge_length + c.edge_length + c.min_below[0]
                if dist < node.min_above[0]:
                    node.min_above = (dist,c.min_below[1])
        # min distance through grandparent
        if not c.parent.is_root():
            dist = node.edge_length + node.parent.min_above[0]
            if dist < node.min_above[0]:
                node.min_above = (dist,node.parent.min_above[1])

    # set up Disjoint Set
    ds = DisjointSet(leaves)
    for node in tree.traverse_preorder(leaves=False):
        # children to min above
        for c in node.children:
            if c.min_below[0] + c.edge_length + node.min_above[0] <= threshold:
                ds.union(c.min_below[1], node.min_above[1])
        for i in range(len(node.children)-1):
            c1 = node.children[i]
            for j in range(i+1, len(node.children)):
                c2 = node.children[j]
                if c1.min_below[0] + c1.edge_length + c2.min_below[0] + c2.edge_length <= threshold:
                    ds.union(c1.min_below[1], c2.min_below[1])
    return [list(s) for s in ds.sets()]

# min_clusters_threshold_max, but all clusters must define a clade
def min_clusters_threshold_max_clade(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()
    for node in tree.traverse_postorder():
        # if I've already been handled, ignore me
        if node.DELETED:
            continue

        # find my undeleted max distances to leaf
        if node.is_leaf():
            node.left_dist = 0; node.right_dist = 0
        else:
            children = list(node.children)
            if children[0].DELETED and children[1].DELETED:
                cut(node); continue
            if children[0].DELETED:
                node.left_dist = 0
            else:
                node.left_dist = max(children[0].left_dist,children[0].right_dist) + children[0].edge_length
            if children[1].DELETED:
                node.right_dist = 0
            else:
                node.right_dist = max(children[1].left_dist,children[1].right_dist) + children[1].edge_length

            # if my kids are screwing things up, cut both
            if node.left_dist + node.right_dist > threshold:
                cluster_l = cut(children[0])
                node.left_dist = 0
                cluster_r = cut(children[1])
                node.right_dist = 0

                # add cluster
                for cluster in (cluster_l,cluster_r):
                    if len(cluster) != 0:
                        clusters.append(cluster)
                        for leaf in cluster:
                            leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# pick the threshold between 0 and "threshold" that maximizes number of (non-singleton) clusters
def argmax_clusters(method,tree,threshold,support):
    from copy import deepcopy
    assert threshold > 0, "Threshold must be positive"
    #thresholds = pairwise_dists_below_thresh(deepcopy(tree),threshold)
    thresholds = [i*threshold/NUM_THRESH for i in range(NUM_THRESH+1)]
    best = None; best_num = -1; best_t = -1
    for i,t in enumerate(thresholds):
        print("%s%%"%str(i*100/len(thresholds)).rstrip('0'),end='\r',file=stderr)
        clusters = method(deepcopy(tree),t,support)
        num_non_singleton = len([c for c in clusters if len(c) > 1])
        if num_non_singleton > best_num:
            best = clusters; best_num = num_non_singleton; best_t = t
    print("\nBest Threshold: %f"%best_t,file=stderr)
    return best

# cut all branches longer than the threshold
def length(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()
    for node in tree.traverse_postorder():
        # if I've already been handled, ignore me
        if node.DELETED:
            continue

        # if i'm screwing things up, cut me
        if node.edge_length is not None and node.edge_length > threshold:
            cluster = cut(node)
            if len(cluster) != 0:
                clusters.append(cluster)
                for leaf in cluster:
                    leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# same as length, and clusters must define a clade
def length_clade(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()
    for node in tree.traverse_postorder():
        # if I've already been handled, ignore me
        if node.DELETED or node.is_leaf():
            continue

        # if either kid is screwing things up, cut both
        children = list(node.children)
        if children[0].edge_length > threshold or children[1].edge_length > threshold:
            cluster_l = cut(children[0])
            cluster_r = cut(children[1])

            # add clusters
            for cluster in (cluster_l,cluster_r):
                if len(cluster) != 0:
                    clusters.append(cluster)
                    for leaf in cluster:
                        leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# cut tree at threshold distance from root (clusters will be clades by definition) (ignores support threshold if branch is below cutting point)
def root_dist(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()
    for node in tree.traverse_preorder():
        # if I've already been handled, ignore me
        if node.DELETED:
            continue
        if node.is_root():
            node.root_dist = 0
        else:
            node.root_dist = node.parent.root_dist + node.edge_length
        if node.root_dist > threshold:
            cluster = cut(node)
            if len(cluster) != 0:
                clusters.append(cluster)
                for leaf in cluster:
                    leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# cut tree at threshold distance from the leaves (if tree not ultrametric, max = distance from furthest leaf from root, min = distance from closest leaf to root, avg = average of all leaves)
def leaf_dist(tree,threshold,support,mode):
    modes = {'max':max,'min':min,'avg':avg}
    assert mode in modes, "Invalid mode. Must be one of: %s" % ', '.join(sorted(modes.keys()))
    dist_from_root = modes[mode](d for u,d in tree.distances_from_root(internal=False)) - threshold
    return root_dist(tree,dist_from_root,support)
def leaf_dist_max(tree,threshold,support):
    return leaf_dist(tree,threshold,support,'max')
def leaf_dist_min(tree,threshold,support):
    return leaf_dist(tree,threshold,support,'min')
def leaf_dist_avg(tree,threshold,support):
    return leaf_dist(tree,threshold,support,'avg')

METHODS = {
    'max': min_clusters_threshold_max,
    'max_clade': min_clusters_threshold_max_clade,
    'sum_branch': min_clusters_threshold_sum_bl,
    'sum_branch_clade': min_clusters_threshold_sum_bl_clade,
    'avg_clade': min_clusters_threshold_avg_clade,
    'med_clade': min_clusters_threshold_med_clade,
    'single_linkage': single_linkage_cut,
    'single_linkage_cut': single_linkage_cut,
    'single_linkage_union': single_linkage_union,
    'length': length,
    'length_clade': length_clade,
    'root_dist': root_dist,
    'leaf_dist_max': leaf_dist_max,
    'leaf_dist_min': leaf_dist_min,
    'leaf_dist_avg': leaf_dist_avg
}
THRESHOLDFREE = {'argmax_clusters':argmax_clusters}

def run_TreeCluster(threshold, tree_file, threshold_free, method, support):
    trees = []
    trees.append(read_tree_newick(tree_file))
    # run algorithm
    for t,tree in enumerate(trees):
        if threshold_free is None:
            clusters = METHODS[method.lower()](tree,threshold,support)
        else:
            clusters = THRESHOLDFREE[threshold_free](METHODS[method.lower()],tree,threshold,support)
        
    return clusters


class StopGridSearch( Exception ):
    pass

def precal_GridSearch(threshold, tree_file, threshold_free, method, support, expected_num_clusters):
        cluster = run_TreeCluster(threshold, tree_file, threshold_free, method, support)
        loss = expected_num_clusters-len(cluster)
        grad = 0.1 * pow(10,-len(str(loss)))
        min_loss= expected_num_clusters

        return loss,grad,min_loss

def GridSearch(loss, grad, min_loss,threshold, tree_file, threshold_free, method, support, max_iter, expected_num_clusters):
        threshold_vals= [threshold]
        loss_vals = [loss]
        current_grad = [grad]
        i=1
        best_val = threshold
        curr_loss= loss
        #print("Iteration : ", 0, " Loss: ", loss, " Threshold_val:", threshold)
        try:
            while(i<=max_iter+1):
                curr_grad = ( 0.1 * pow(10,-len(str(loss_vals[i-1]))) ) 
                #print(curr_grad)
                curr_val = threshold_vals[i-1] - (grad*curr_loss)*np.sign(grad*loss)
                if curr_val<0:
                    raise StopGridSearch
                    curr_val = best_val
                    min_loss = loss_vals[i]
                curr_loss = abs(expected_num_clusters-len(run_TreeCluster(curr_val, tree_file, threshold_free, method, support)))
                threshold_vals.append(curr_val)
                loss_vals.append(curr_loss)
                #print(loss_vals[i], min_loss)
                if loss_vals[i] < min_loss:
                    min_loss = loss_vals[i] 
                    best_val = curr_val
        
                #print(" Iteration : ", i, " Loss: ", curr_loss, " Threshold_val:", curr_val, "Min_loss:",min_loss )
         
                i+=1
        
        except StopGridSearch:
            pass
        
        return best_val, min_loss
    


def someiter_search(tree_file,threshold_val,threshold_free,method,support,expected_num_clusters,max_iter,max_iter_cont):
   j=0
   while j<=max_iter_cont:
        if j==0:
            loss,grad,min_loss = precal_GridSearch(threshold_val, tree_file, threshold_free, method, support, expected_num_clusters)
            best_min_loss = min_loss
            curr_threshold_val = threshold_val
            best_threshold_val = threshold_val
        else:
            #print(loss, grad, min_loss,curr_threshold_val, tree_file, threshold_free, method, support, max_iter, expected_num_clusters)
            curr_threshold_val, min_loss = GridSearch(loss, grad, min_loss,curr_threshold_val, tree_file, threshold_free, method, support, max_iter, expected_num_clusters)
            grad = min_loss * ( 0.1 * pow(10,-len(str(min_loss))) ) * pow(.1,j)
            loss = min_loss
            if loss< best_min_loss:
                best_min_loss= min_loss
                best_threshold_val = curr_threshold_val

        #print("ITER:",j, " Loss: ", best_min_loss, " Threshold_val:", best_threshold_val)
        j+=1

   return best_min_loss, best_threshold_val 


def search_best_thresh(threshold_val_iter,max_iter,max_iter_cont,expected_num_clusters,tree_file,threshold_free,method,support):
#threshold_val_iter,expected_num_clusters,tree_file="phylogeny.nwk",threshold_free=None,method='avg_clade',support=float('-inf'),max_iter = 10):
    thresh_list = list(np.random.uniform(0,1,threshold_val_iter))

    for val in thresh_list:
        best_min_loss = expected_num_clusters
        
        min_loss, curr_threshold_val = someiter_search(tree_file,val,threshold_free,method,support,expected_num_clusters,max_iter,max_iter_cont)
        if min_loss < best_min_loss:
                best_min_loss= min_loss
                best_threshold_val = curr_threshold_val
            
    
        print("REAL ITER:",thresh_list.index(val)+1, " Loss: ", min_loss, " Threshold_val:", curr_threshold_val)
    
    return best_threshold_val



def return_best_cluster(initvalue,coarsevalue,finevalue,n_clust,inp,tf,method,support):

    val = search_best_thresh(initvalue,coarsevalue,finevalue,n_clust,inp,tf,method,support)

    return run_TreeCluster(val, inp, tf,  method, support), val

if __name__ == "__main__":
    # parse user arguments
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input Tree File")
    parser.add_argument('-s', '--support', required=False, type=float, default=float('-inf'), help="Branch Support Threshold")
    parser.add_argument('-it', '--initvalue', required=False, type=int, default=int(10), help="Threshold Value Initalization for Regressions")
    parser.add_argument('-co', '--coarsevalue', required=False, type=int, default=int(40), help="Coarse Value adjustment for Regressions")
    parser.add_argument('-n', '--n_clust', required=True, type=int,  help="Number of Clusters")
    parser.add_argument('-fi', '--finevalue', required=False, type=int, default=int(10), help="Fine Value adjustment for Regressions")
    parser.add_argument('-m', '--method', required=False, type=str, default='max_clade', help="Clustering Method (options: %s)" % ', '.join(sorted(METHODS.keys())))
    #parser.add_argument('-tf', '--threshold_free', required=False, type=str, default=None, help="Threshold-Free Approach (options: %s)" % ', '.join(sorted(THRESHOLDFREE.keys())))
    args = parser.parse_args()
    assert args.method.lower() in METHODS, "ERROR: Invalid method: %s" % args.method
    #assert args.threshold_free is None or args.threshold_free in THRESHOLDFREE, "ERROR: Invalid threshold-free approach: %s" % args.threshold_free
    assert args.support >= 0 or args.support == float('-inf'), "ERROR: Branch support must be at least 0"
    assert args.initvalue >=0, "ERROR: Number of Theeshold Value Init Regressions must at least be 1"
    assert args.n_clust >=0, "ERROR: Number of Clusters must at least be 1"
    assert args.finevalue >=0, "ERROR: Number of  Value adjustment Regressions must at least be 1"
    assert args.coarsevalue >=0, "ERROR: Number of  Value adjustment Regressions must at least be 1"

    cluster = return_best_cluster(args.initvalue,args.coarsevalue,args.finevalue,args.n_clust,args.input,None,args.method,args.support)


    print(cluster)
