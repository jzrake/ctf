#!/usr/bin/env python


def CopyH5File(src, dst, skip):
    from subprocess import Popen, PIPE
    from os import system
    p = Popen("h5ls %s" % src, shell=True, stdout=PIPE, stderr=PIPE)
    output, errors = p.communicate()
    contents = [line.split()[0] for line in output.split("\n") if line]
    system("rm -f %s" % dst)

    for entry in contents:
        if entry in skip: continue

        params = {"Infile": src,
                  "Outfile": dst,
                  "ObjectName": entry}
        cmd = \
            "h5copy -v -i %(Infile)s -o %(Outfile)s " \
            "-s %(ObjectName)s -d %(ObjectName)s" % params

        system(cmd)


if __name__ == "__main__":
    from optparse import OptionParser
    from h5py import File
    from numpy import array, zeros
    from itertools import product

    parser = OptionParser()
    opts, args = parser.parse_args()

    CopyH5File(src=args[0], dst=args[1], skip='prim')

    h5filei = File(args[0], 'r')
    h5fileo = File(args[1], 'a')

    new_dt =  h5fileo['status/Timestep'].value / 2.0
    del h5fileo['status/Timestep']
    h5fileo['status/Timestep'] = new_dt

    primi = h5filei["prim"]
    primo = h5fileo.create_group("prim")

    for dset in primi:
        print "Working on data set", dset
        datai = primi[dset]

        if not datai.chunks:
            chnksize = array([16,16,16])
        else:
            chnksize = array(datai.chunks)

        globsize = array(datai.shape)
        numchunk = globsize/chnksize
        datao = primo.create_dataset(dset, globsize*2, chunks=tuple(chnksize*2), dtype='f8')
        print "Using chunk size", datao.chunks

        # Below, 'ijk' labels the index of the chunk. IJK labels the absolute
        # index of the first zone in that chunk, for the input array. That data
        # is buffered in memory for re-writing into the upsampled array. The
        # nested loop over the iterator 'bit' takes care of moving over the 8
        # zones in the finer array which will receive the value of the single
        # zone in the courser array.
        for ijk in product(*(range(n) for n in numchunk)):
            print "Loading chunk", ijk
            IJK = ijk * chnksize
            buf = datai[IJK[0]:(IJK+chnksize)[0],
                        IJK[1]:(IJK+chnksize)[1],
                        IJK[2]:(IJK+chnksize)[2]]
            buf8 = zeros([s*2 for s in buf.shape])

            for bit in product(*(range(2),)*3):
                buf8[bit[0]::2, bit[1]::2, bit[2]::2] = buf

            datao[IJK[0]*2:(IJK+chnksize)[0]*2,
                  IJK[1]*2:(IJK+chnksize)[1]*2,
                  IJK[2]*2:(IJK+chnksize)[2]*2] = buf8
