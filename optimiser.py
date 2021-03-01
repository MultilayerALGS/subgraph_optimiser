# -*- coding: utf-8 -*-

import logging
from math import log as ln
import numpy as np
import copy
from utils import EmptyGraph, powerset

formatter = logging.Formatter('%(asctime)s %(message)s', '%H:%M:%S')
LOG = logging.getLogger(__name__)
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(formatter)
LOG.addHandler(consoleHandler)
#LOG.setLevel(logging.INFO)



class Optimiser:
    def __init__(self, adj_matrix, data, plus_fun="max", minus_fun="max", fails_allowed=0):
        self.graph = EmptyGraph(adj_matrix)
        self.data = data
        self.plus_fun = plus_fun
        self.minus_fun = minus_fun
        self.fails_allowed = fails_allowed

    def best_possible(self):
        smallest_diff = None
        for u in self.graph.vertices():
            for v in self.graph.possNbs(u):
                diff  = abs(self.data[u] - self.data[v])
                if smallest_diff is None or diff < smallest_diff:
                    smallest_diff = diff
        return (1/2 * sum(ln(self.graph.possibleDegree(u)) for u in self.graph.vertices())
                - (self.graph.size()/2) * ln(2*smallest_diff**2))

    def cond_sorted_nbs_remove(self, i,nbs,options,avs,n):
        fixed_nbs = [nb for nb in nbs if nb not in options]
        nbr_options = [fixed_nbs + poss_nbs for poss_nbs in powerset(options) if len(fixed_nbs + poss_nbs) > 0]
        nbr_list = [(nbs,self.new_vertex_val(i,nbs,avs,n)) for nbs in nbr_options]
        return sorted(nbr_list,key=lambda x: -x[1])

    def cond_sorted_nbs_add(self, i,nbs,options,avs,n):
        fixed_nbs = list(self.graph.nbs(i))
        nbr_options = [fixed_nbs + poss_nbs for poss_nbs in powerset(options) if len(fixed_nbs + poss_nbs) > 0]
        nbr_list = [(nbs,self.new_vertex_val(i,nbs,avs,n)) for nbs in nbr_options]
        return sorted(nbr_list,key=lambda x: -x[1])

    def cond_sorted_nbs_mult(self, i,nbs,options,data,avs,n):
        fixed_nbs = [nb for nb in nbs if nb not in options]
        nbr_options = [fixed_nbs + poss_nbs for poss_nbs in powerset(options) if len(fixed_nbs + poss_nbs) > 0]
        nbr_list = [(nbs,vertex_val_multiple(i,nbs,data,avs,n)) for nbs in nbr_options]
        return sorted(nbr_list,key=lambda x: -x[1])

    def find_best_cond(self, nbr_list,with_nbs,without_nbs):
        i = 0
        found = False
        while (i < len(nbr_list)):
            cand_nbs = nbr_list[i][0]
            found = True
            for w in with_nbs:
                if w not in cand_nbs:
                    found = False
            for w in without_nbs:
                if w in cand_nbs:
                    found = False
            if found:
                return nbr_list[i]
            i = i + 1
        print(f"{with_nbs=} {without_nbs=} {nbr_list=}")
        raise Exception("No best found")

    def find_worst_cond(self, nbr_list,with_nbs,without_nbs):
        found = False
        for cand_pair in reversed(nbr_list):
            cand_nbs = cand_pair[0]
            found = True
            for w in with_nbs:
                if w not in cand_nbs:
                    found = False
            for w in without_nbs:
                if w in cand_nbs:
                    found = False
            if found:
                return cand_pair

    def find_avg_cond(self, nbr_list,with_nbs,without_nbs):
        found = False
        total = 0
        count = 0
        for cand_pair in reversed(nbr_list):
            cand_nbs = cand_pair[0]
            found = True
            for w in with_nbs:
                if w not in cand_nbs:
                    found = False
            for w in without_nbs:
                if w in cand_nbs:
                    found = False
            if found:
                total += cand_pair[1]
                count += 1
        return total / count

    def greedy_opt_remove(self, avsum, order, rev_order, only_later_edges, must_not_remove):
        n = self.graph.size()
        removed = []
        # for every vertex in the matrix
        for vx_ind in range(0, n):
            vx = order[vx_ind]
            try:
                this_nbs =  self.graph.nbs(vx)
            except IndexError as e:
                print(f"vx is {vx}")
            LOG.debug(f"On {vx=} {self.graph.degree(vx)=}")
            if self.graph.degree(vx) <= 1:
                continue
            if only_later_edges:
                nb_options = [nb for nb in self.graph.nbs(vx) if rev_order[nb] > vx_ind and self.graph.degree(nb) > 1 and [vx, nb] not in must_not_remove]
            else:
                nb_options = [nb for nb in self.graph.nbs(vx) if self.graph.degree(nb) > 1 and [vx, nb] not in must_not_remove and [nb, vx] not in must_not_remove]
            # nb_options is the set of neighbours from which we can potentially delete - those that have degree greater than one,
            # and which we haven't yet considered
            sorted_nbd_list = self.cond_sorted_nbs_remove(vx, self.graph.nbs(vx), nb_options,avsum,n)
            # sorted_nbd_list is the sorted list of all neighbour options for vx
            LOG.debug(f"{sorted_nbd_list=}")
            nbr_sorted_nbd_lists = []
            for u in nb_options:
                if only_later_edges:
                    u_options = [w for w in self.graph.nbs(u) if rev_order[w] >= vx_ind and self.graph.degree(w) > 1 and [u, w] not in must_not_remove]
                else:
                    u_options = [w for w in self.graph.nbs(u) if self.graph.degree(w) > 1 and [u, w] not in must_not_remove and [w, u] not in must_not_remove]
                # u_options is the set of neighbours for u that we could still change
                nbr_sorted_nbd_lists += [self.cond_sorted_nbs_remove(u,self.graph.nbs(u),u_options,avsum,n)]
                # add the sorted list of neighbour options for u to the nbr_sorted_nbd_lists
            LOG.info(f"{vx=} nbrs={self.graph.degree(vx)} nb_options={len(nb_options)} sort_list={len(sorted_nbd_list)} sum_nbrs={sum(len(x) for x in nbr_sorted_nbd_lists)}")
            remove = []
            for i in range(0,len(nb_options)):
                nb = nb_options[i]
                LOG.debug(f"\t{nb=}")
                if self.plus_fun == "max":
                    pos_fun = lambda x,y,z: self.find_best_cond(x,y,z)[1]
                elif self.plus_fun == "min":
                    pos_fun = lambda x,y,z: self.find_worst_cond(x,y,z)[1]
                elif self.plus_fun == "avg":
                    pos_fun = lambda x,y,z: self.find_avg_cond(x,y,z)
                if self.minus_fun == "max":
                    neg_fun = lambda x,y,z: self.find_best_cond(x,y,z)[1]
                elif self.minus_fun == "min":
                    neg_fun = lambda x,y,z: self.find_worst_cond(x,y,z)[1]
                elif self.minus_fun == "avg":
                    neg_fun = lambda x,y,z: self.find_avg_cond(x,y,z)
                vx_gain = pos_fun(sorted_nbd_list,[nb],[]) - neg_fun(sorted_nbd_list,[],[nb])
                LOG.debug(f"\tvx gain: {pos_fun(sorted_nbd_list, [nb],[])} - {neg_fun(sorted_nbd_list, [],[nb])}")
                LOG.debug(f"\tnb gain: {pos_fun(nbr_sorted_nbd_lists[i], [vx],[])} - {neg_fun(nbr_sorted_nbd_lists[i], [],[vx])}")
                nb_gain = pos_fun(nbr_sorted_nbd_lists[i],[vx],[]) - neg_fun(nbr_sorted_nbd_lists[i],[],[vx])
                combined_gain = vx_gain + nb_gain
                LOG.debug(f"\ttry {nb=}, {combined_gain=}, {vx_gain=}, {nb_gain=}")
                #LOG.debug(f"\ttry keeping {find_best_cond(sorted_nbd_list,[],[nb])[0]},{vx} and {find_best_cond(nbr_sorted_nbd_lists[i],[],[vx])[0]},{nb}")
                if (combined_gain < 0):
                    remove.append([[vx,nb],combined_gain])
            LOG.debug(f"\t\t{remove=}")
            if remove:
                self.graph.remove_edges([e[0] for e in remove])
                removed.extend(e[0] for e in remove)
        return removed

    def greedy_opt_add(self, avsum, order, rev_order, only_later_edges, must_not_add):
        n = self.graph.size()
        added = []
        # for every vertex in the matrix
        for vx_ind in range(0, n):
            vx = order[vx_ind]
            LOG.debug(f"On {vx=}")
            if only_later_edges:
                nb_options = [nb for nb in self.graph.possNbs(vx) if rev_order[nb] > vx_ind and [vx, nb] not in must_not_add]
            else:
                nb_options = [nb for nb in self.graph.possNbs(vx) if [vx, nb] not in must_not_add and [nb, vx] not in must_not_add]
            # nb_options is the set of neighbours from which we can potentially delete - those that have degree greater than one,
            # and which we haven't yet considered
            sorted_nbd_list = self.cond_sorted_nbs_add(vx, self.graph.allPossNbs(vx), nb_options,avsum,n)
            # sorted_nbd_list is the sorted list of all neighbour options for vx
            LOG.debug(f"{sorted_nbd_list=}")
            nbr_sorted_nbd_lists = []
            for u in nb_options:
                if only_later_edges:
                    u_options = [w for w in self.graph.possNbs(u) if rev_order[w] >= vx_ind]
                else:
                    u_options = self.graph.possNbs(u)
                # u_options is the set of neighbours for u that we could still change
                nbr_sorted_nbd_lists += [self.cond_sorted_nbs_add(u,self.graph.allPossNbs(u),u_options,avsum,n)]
                # add the sorted list of neighbour options for u to the nbr_sorted_nbd_lists
            add = []
            LOG.info(f"{vx=} nbrs={self.graph.degree(vx)} nb_options={len(nb_options)} sort_list={len(sorted_nbd_list)} sum_nbrs={sum(len(x) for x in nbr_sorted_nbd_lists)}")
            for i in range(0,len(nb_options)):
                nb = nb_options[i]
                LOG.debug(f"\t{nb=}")
                if self.plus_fun == "max":
                    pos_fun = lambda x,y,z: self.find_best_cond(x,y,z)[1]
                elif self.plus_fun == "min":
                    pos_fun = lambda x,y,z: self.find_worst_cond(x,y,z)[1]
                elif self.plus_fun == "avg":
                    pos_fun = lambda x,y,z: self.find_avg_cond(x,y,z)
                if self.minus_fun == "max":
                    neg_fun = lambda x,y,z: self.find_best_cond(x,y,z)[1]
                elif self.minus_fun == "min":
                    neg_fun = lambda x,y,z: self.find_worst_cond(x,y,z)[1]
                elif self.minus_fun == "avg":
                    neg_fun = lambda x,y,z: self.find_avg_cond(x,y,z)
                vx_gain = pos_fun(sorted_nbd_list,[nb],[]) - neg_fun(sorted_nbd_list,[],[nb])
                LOG.debug(f"\tvx gain: {pos_fun(sorted_nbd_list, [nb],[])} - {neg_fun(sorted_nbd_list, [],[nb])}")
                LOG.debug(f"\tnb gain: {pos_fun(nbr_sorted_nbd_lists[i], [vx],[])} - {neg_fun(nbr_sorted_nbd_lists[i], [],[vx])}")
                nb_gain = pos_fun(nbr_sorted_nbd_lists[i],[vx],[]) - neg_fun(nbr_sorted_nbd_lists[i],[],[vx])
                combined_gain = vx_gain + nb_gain
                LOG.debug(f"\ttry {nb=}, {combined_gain=}, {vx_gain=}, {nb_gain=}")
                #LOG.debug(f"\ttry keeping {find_best_cond(sorted_nbd_list,[],[nb])[0]},{vx} and {find_best_cond(nbr_sorted_nbd_lists[i],[],[vx])[0]},{nb}")
                if (combined_gain > 0):
                    add.append([[vx,nb],combined_gain])
            LOG.debug(f"\t\t{add=}")
            if add:
                self.graph.add_edges([e[0] for e in add])
                added.extend(e[0] for e in add)
        return added



    def greedy_opt_multiple(matrix,data,avsums):
        n = len(matrix)
        newmatrix = matrix
        for vx_ind in range(0, len(matrix)):
            vx = order[vx_ind]
            if len(get_nbs(vx,newmatrix)) > 1:
                nb_options = [nb for nb in get_nbs(vx,newmatrix) if rev_order(nb) > vx_ind and len(get_nbs(nb,newmatrix)) > 1]
                sorted_nbd_list = cond_sorted_nbs_mult(vx,get_nbs(vx,newmatrix),nb_options,data,avsums,n)
                nbr_sorted_nbd_lists = []
                for u in nb_options:
                    u_options = [w for w in get_nbs(u,newmatrix) if w >= vx and len(get_nbs(w,newmatrix)) > 1]
                    nbr_sorted_nbd_lists += [cond_sorted_nbs_mult(u,get_nbs(u,newmatrix),u_options,data,avsums,n)]
                remove = []
                for i in range(0,len(nb_options)):
                    nb = nb_options[i]
                    vx_gain = find_best_cond(sorted_nbd_list,[nb],[])[1] - find_best_cond(sorted_nbd_list,[],[nb])[1]
                    nb_gain = find_best_cond(nbr_sorted_nbd_lists[i],[vx],[])[1] - find_best_cond(nbr_sorted_nbd_lists[i],[],[vx])[1]
                    combined_gain = vx_gain + nb_gain
                    if (combined_gain < 0):
                        remove.append([[vx,nb],combined_gain])
                if len(remove) < len(get_nbs(vx,newmatrix)):
                    newmatrix = remove_edges(newmatrix,[e[0] for e in remove])
                else:
                    sorted_remove = sorted(remove,key=lambda x: x[1])
                    newmatrix = remove_edges(newmatrix,[e[0] for e in sorted_remove[1:]])
        return newmatrix


    def default_avsum(self):
        return sum([self.disc(i) for i in range(0,self.graph.size())])


    def new_vertex_val(self, index,nbs,avsum,n):
        # deg(index) * ND(index)
        vx_av = len(nbs)*(self.data[index] - sum([self.data[u] for u in nbs])/len(nbs))**2
        if vx_av > avsum:
            # the discrepancy in the proposed graph at this vertex alone is already
            # worse than the total discrepancy in the old. As long as we are only
            # removing edges, this increase in penalty cannot be counter-acted, so
            # skip some invalid logs of negatives, and just return a really big
            # negative number
            return -1e20
        ub = np.log(len(nbs))/2 - (n/2)*(np.log1p(vx_av/(avsum - vx_av)))
        return ub


    def vertex_val_multiple(self, index,nbs,avsums,n):
        contrib = 0
        deg = len(nbs)
        q = len(self.data)
        for i in range(0,q):
            vx_av = len(nbs)*(self.data[i][index] - sum([self.data[i][u] for u in nbs])/deg)**2
            val = q*np.log(deg)/2 - (n*q/2)*(np.log(1 + vx_av/max(0.001,avsums[i]-vx_av)))
            contrib += val
        return contrib


    def disc(self, index):
        deg = self.graph.degree(index)
        result = deg*(self.data[index] - sum([self.data[j] for j in self.graph.nbs(index)])/deg)**2
        return result

    def matrix_score(self):
        n = self.graph.size()
        pluspart = sum([np.log(self.graph.degree(i)) for i in range(0,n)])/2
        penalty = (n/2)*np.log(sum([self.disc(index) for index in range(0,n)])/n)
        return pluspart - penalty


    def matrix_score_mult(self):
        n = self.graph.size()
        q = len(self.data)
        pluspart = (q/2)*sum([np.log(self.graph.degree(i)) for i in range(0,n)])
        penalty = (n*q/2)*np.log(sum([sum([self.disc(index) for index in range(0,n)]) for values in self.data])/(n*q))
        return pluspart- penalty

    def iterative_opt(self, order=None, rev_order=None, only_later_edges=True, starter_edges=None,
                      remove=True, add=True, remove_first=False):
        all_added = []
        if order is None:
            order = list(x for x in range(len(self.data)))
            rev_order = list(x for x in range(len(self.data)))
        if not add or (remove and remove_first):
            # If not adding , then we start with full graph
            self.graph.add_edges(self.graph.origEdges())
            LOG.info("Removing from original graph")
        else:
            LOG.info("Building up from empty graph")
            # We are building from a mostly empty graph
            if starter_edges:
                for edge in starter_edges:
                    self.graph.add_edge(edge)
            for u in range(self.graph.size()):
                if self.graph.degree(u) >= 1:
                    continue
                best = 99999999
                bestv = None
                for v in self.graph.possNbs(u):
                    if u == v:
                        continue
                    diff = abs(self.data[u] - self.data[v])
                    if bestv is None or diff < best:
                        bestv = v
                        best = diff
                self.graph.add_edge([u,bestv])
        must_not_add = []
        must_not_remove = []
        if starter_edges:
            must_not_remove.extend(starter_edges)
            print(f"{must_not_remove=}")
        oldscore = self.matrix_score()
        LOG.info(f"Run starting:{add=} {remove=} {remove_first=} {oldscore=:.3f}")
        change = True
        while change:
            change = False
            if add and not (remove_first and remove):
                mode = "add"
                def_avsum = self.default_avsum()
                oldscore = self.matrix_score()
                lastadded = self.greedy_opt_add(def_avsum, order, rev_order, only_later_edges, must_not_add)
                newscore = self.matrix_score()
                if lastadded:
                    if newscore < oldscore:
                        res = "bad"
                        must_not_add.extend(lastadded)
                        self.graph.remove_edges(lastadded)
                    else:
                        res = "good"
                        change = True
                else:
                    res = "---"
                change_size = len(lastadded)
                LOG.info(f"{mode}\t{res}\t{oldscore=:.3f}\t{newscore=:.3f} \t{change_size=}")
            remove_first = False
            if remove:
                mode = "remove"
                def_avsum = self.default_avsum()
                oldscore = self.matrix_score()
                lastremoved = self.greedy_opt_remove(def_avsum, order, rev_order, only_later_edges, must_not_remove)
                newscore = self.matrix_score()
                if lastremoved:
                    if newscore < oldscore:
                        res = "bad"
                        must_not_remove.extend(lastremoved)
                        self.graph.add_edges(lastremoved)
                    else:
                        change = True
                        res = "good"
                else:
                    res = "---"
                change_size = len(lastremoved)
                LOG.info(f"{mode}\t{res}\t{oldscore=:.3f}\t{newscore=:.3f} \t{change_size=}")
        return self.graph


    def iterative_opt_multiple(matrix,data):
        newscore = matrix_score_mult(matrix,data)
        oldscore = -100000
        while oldscore < newscore:
            oldscore = newscore
            def_avsums = [default_avsum(oldmatrix,dat) for dat in data]
            newmatrix = self.greedy_opt_multiple(def_avsums)
            newscore = self.matrix_score_mult()
        raise Exception("This is broken, need to re-add removed edges from last iteration")
        return oldmatrix


