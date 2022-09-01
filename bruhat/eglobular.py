#!/usr/bin/env python

from typing import NamedTuple, Dict,List, Optional

# AAARRRGGGGGGGGHHHH



#class ENode(NamedTuple):
#    name: str
#    args: 'Tuple[EClassID, ...]'

class ENode(object):

    def __init__(self, name, args=()):
        self.name = str(name)
        for arg in args:
            assert isinstance(arg, EClassID), repr(arg)
        self.args = tuple(args)

    def __eq__(self, other):
        return (self.name, self.args) == (other.name, other.args)

    def __hash__(self):
        return hash((self.name, self.args))

    def __mul__(self, other):
        return ENode("*", (self, other))

    def canonicalize(self): 
        return ENode(self.name, tuple(arg.find() for arg in self.args))
        

class EClassID(object):
    def __init__(self, egraph, id):
        self.id = id #  just for debugging
        self.egraph = egraph # not used ??
        #self.parent : 'Optional[EClassID]' = None
        self.parent = None # EClassID

        # A list of enodes that use this eclass and the EClassID of that use.
        # This is only set on a canonical EClassID (parent == None) and
        # will be used to repair the graph on merges.
        #self.uses : Optional[List[Tuple['ENode', 'EClassID']]] = []
        self.uses = [] # list of (ENode, EClassID)

    def __repr__(self):
        return f'e{self.id}'

    # union-find's 'find' operation, that finds the canonical object for this set
    def find(self):
        if self.parent is None:
            return self
        r = self.parent.find()
        # path compression, makes find cheaper on subsequent calls.
        self.parent = r
        return r

    def __hash__(self):
        return id(self)

    def __eq__(lhs, rhs):
        lhs = lhs.find()
        rhs = rhs.find()
        return lhs is rhs

Env = Dict[str, EClassID]

class EGraph(object):
    def __init__(self):
        self.i = 0 # id counter for debugging EClassID's
        self.version = 0 # increments every time we mutate the EGraph, so that we
                         # can quickly see when we haven't changed one

        # this map, from ENodes (canonical) to EClassID (not canonical), is maintained
        # so that we check if we have already defined a particular ENode
        self.hashcons : Dict['ENode', 'EClassID'] = {}

        # a list of EClasseID's that have mutated by a merge operation, and whose users must
        # be updated to a new canonical version, we will later define a `rebuild` function
        # that processes this worklist.
        self.worklist : List['EClassID'] = []


    def add_enode(self, enode: 'ENode'):
        enode = enode.canonicalize()
        eclass_id = self.hashcons.get(enode, None)
        if eclass_id is None:
            # we don't already have this node, so create one
            self.version += 1
            eclass_id = self._new_singleton_eclass()
            for arg in enode.args:
                # we need to update the uses lists of each arg,
                # since we are adding this enode as a new use.
                arg.uses.append((enode, eclass_id))
            self.hashcons[enode] = eclass_id
        # hashcons's rhs is not canonicalized, so canonicalize it here:
        return eclass_id.find()

    def _new_singleton_eclass(self):
        r = EClassID(self, self.i)
        self.i += 1
        return r

    # extract the actual EClassID's
    def eclasses(self) -> Dict['EClassID', List['ENode']]:
        r = {}
        for node, eid in self.hashcons.items():
            eid = eid.find()
            if eid not in r:
                r[eid] = [node]
            else:
                r[eid].append(node)
        return r

    #    Merging EClassID'es
    #    ================
    #    
    #    Notice that `e2` and `e5` can be merged, since they define
    #    the same value. To do this we need to define a the `merge`
    #    operation. The details are described in the paper, but
    #    you can follow the rest of this document without understanding all the details.
    
    # tell the graph that 'a' and 'b' calculate the same value.
    def merge(self, a: 'EClassID', b: 'EClassID'):
        a = a.find()
        b = b.find()
        if a is b:
            return a
        self.version += 1
        a.parent = b
        # maintain the invariant that uses are only
        # recorded on the top level EClassID
        b.uses += a.uses
        a.uses = None
        # we have updated eclass b, so nodes in the hashcons
        # might no longer be canonical and we might discover that two
        # enodes are actually the same value. We will repair this later, by
        # remember what EClassID's changed:
        self.worklist.append(b)
    
    # ensure we have a de-duplicated version of the EGraph
    def rebuild(self):
        while self.worklist:
            # deduplicate repeated calls to repair the same eclass
            todo = { eid.find(): None for eid in self.worklist}
            self.worklist = []
            for eclass_id in todo:
                self.repair(eclass_id)
    
    def repair(self, eclass_id: 'ENodeID'):
        assert eclass_id.parent is None
        # reset the uses of this eclass, we will fill them in again at the end
        uses, eclass_id.uses = eclass_id.uses, []
        # any of the uses in the hashcons might no longer be canonical, so re-canonicalize it
        for p_node, p_eclass in uses:
            if p_node in self.hashcons:
                del self.hashcons[p_node]
            p_node = p_node.canonicalize()
            self.hashcons[p_node] = p_eclass.find()
        # because we merged classes, some of the enodes that are uses might now be the same expression,
        # meaning we can merge further EClassID's
        new_uses = {}
        for p_node, p_eclass in uses:
            p_node = p_node.canonicalize()
            if p_node in new_uses:
                self.merge(p_eclass, new_uses[p_node])
            new_uses[p_node] = p_eclass.find()
        # note the find, it is possible that eclass_id got merged
        # and uses should only be attributed to the eclass representative
        eclass_id.find().uses += new_uses.items()


def test():

    a = ENode("a")
    aa = a*a
    
    egraph = EGraph()
    enode = egraph.add_enode(aa)

    return
    
    #    We can add other expressions as well, like `a << 1`, which is equivalent to `a * 2`:
    
    egraph.add_node(express(lambda a: a << 1))
    
    #    We can try it out by merging our nodes and observering
    #    that `e5` now has two definitions.
    
    # lookup the EClassID's for these nodes
    x_times_2 = egraph.add_node(express(lambda a: a * 2))
    x_shift_1 = egraph.add_node(express(lambda a: a << 1))
    assert x_times_2 != x_shift_1

    eid = egraph.add_node(express(lambda a: a << 1))
    assert eid is x_shift_1

    # and merge them
    egraph.merge(x_times_2, x_shift_1)
    egraph.rebuild()

    assert x_times_2 == x_shift_1


if __name__ == "__main__":
    test()
    print("OK\n")



