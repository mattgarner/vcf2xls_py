def overlap_alignments(align, overlap):
    print "align", align
    # create a list of starts and ends
    stends = [ (a[0], 1) for a in align ]
    stends += [ (a[1] + overlap, -1) for a in align ]
    stends.sort(key=lambda x: x[0])

    # now we should have a list of starts and ends ordered by position,
    # e.g. if the ranges are 5..10, 8..15, and 12..13, we have
    # (5,1), (8,1), (10,-1), (12,1), (13,-1), (15,-1)

    # next, we form a cumulative sum of this
    s = 0
    cs = []
    print stends
    for se in stends:
        print "se", se
        s += se[1]
        print "s", s
        cs.append((se[0], s))
        print "cs", cs
        print

    print "final cs", cs

    # this is, with the numbers above, (5,1), (8,2), (10,1), (12,2), (13,1), (15,0)
    # so, 5..8 belongs to one range, 8..10 belongs to two overlapping range,
    # 10..12 belongs to one range, etc

    # now we'll find all contiguous ranges
    # when we traverse through the list of depths (number of overlapping ranges), a new
    # range starts when the earlier number of overlapping ranges has been 0
    # a range ends when the new number of overlapping ranges is zero 
    prevdepth = 0
    start = 0
    combined = []
    for pos, depth in cs:
        if prevdepth == 0:
            start = pos
        elif depth == 0:
            combined.append((start, pos-overlap))
        prevdepth = depth

    print combined
    return combined

overlap_alignments( ([5,10],[8,15],[12,13],[16,20]) , 0 )