import sys
import xmlrpc.client as xmlrpclib
import argparse
import pymolr

def load_system(coordx, top, frame, sele, cpmd, dump, is_trr):
    if is_trr:
        if top:
            a = pymolr.System(mpt_file=top, trr_files=coordx)
            print(f"Reading {a.nframes} frames from {coordx}\n")
            if frame == -1: frame = a.nframes - 1
            print(f"You have selected frame no {frame}\n")
            sel = a.select(sele, frame=frame)
            print(f"which is step number {sel.step} from the GROMACS trajectory\n")
            return sel
        else:
            print("Topology not specified!\n")
            sys.exit(1)
    else:
        print(f"Reading single frame from {coordx}\n")
        if not top:
            print("Topology not specified!\n")
            sys.exit(1)
        else:
            return pymolr.select_coords(top, coordx, sele, cpmd, dump)

def process_sele(sele):
    if sele.split('.')[-1] == 'dat':
        with open(sele, 'r') as f:
            return f.read()
    elif sele.split('.')[-1] == 'inp':
        import mimicpy
        cpmd = mimicpy.CpmdScript.from_file(sele)
        ids = [i.split()[1] for i in cpmd.MIMIC.OVERLAPS.splitlines()[1:]]
        return 'id is ' + ' or id is '.join(ids)
    else:
        return sele

def main():
    print('\n \t                ***** PyMOLR *****                  ')
    print('\n \t For more information type pymolr --help \n')
    
    parser = argparse.ArgumentParser(prog='pymolr')
        
    parser.add_argument("coords", help="coordinate file", type=str)
    parser.add_argument("top", help="topology file", type=str, nargs='?', default=None)
    parser.add_argument("-i", "--info", help="print info about trr", default=False, type=bool, metavar='')
    parser.add_argument("-u", "--url", help="PyMOL remote server URL", type=str, metavar='')
    parser.add_argument("-n", "--name", help="name of molecule in PyMOL", metavar='')
    parser.add_argument("-f", "--frame", help="frame to read if trr", metavar='', type=int)
    parser.add_argument("-r", "--ref", help="reference to fix pbc", metavar='', type=str)
    parser.add_argument("-s", "--sele",  nargs='+', help="selection query", metavar='')
    parser.add_argument("-rs", "--rsele",  nargs='+', help="selection query", metavar='')
    parser.add_argument("-c", "--cpmd", help="CPMD input filename", metavar='', type=str)
    parser.add_argument("-d", "--dump", help="TPR dump filename", metavar='', type=str)
    parser.add_argument("-o", "--out", help="output coordinate filename", metavar='', type=str)
    parser.set_defaults(feature=False)
    args = parser.parse_args()
    
    coordx = args.coords
    top = args.top
    frame = 0 if args.frame is None else int(args.frame)
    sele = 'all' if args.sele is None else process_sele(' '.join(args.sele))
    ref = args.ref
    rsele = 'all' if args.rsele is None else process_sele(' '.join(args.rsele))
    out = args.out
    name = 'mol1' if args.name is None else args.name
    info = args.info
    url = 'http://localhost:9123' if args.url is None else args.url
    cpmd = args.cpmd
    dump = args.dump
        
    if coordx.split('.')[-1] == 'trr':
        is_trr = True
    else:
        is_trr = False
        if sele == 'all' and not ref:
            if out:
                print(f"Why do you want to write the contents of {coordx} into {out} without any changes?\nI'm not a copying program!\n")
            pymol = xmlrpclib.ServerProxy(url)
            with open(coordx, 'r') as f:
                s = f.read().replace('\n','@')
            print(f"Connecting to PyMOL at address: {url}")
            pymol.do(f'pymolr_load("{s}", "{name}")')
    
    selected_frame = load_system(coordx, top, frame, sele, cpmd, dump, is_trr)
    
    if ref:
        if not rsele:
            selected_frame.fix_pbc(pymolr.select_coords(top, ref, sele))
        else:
            ref_coords = pymolr.select_coords(top, ref, rsele, cpmd, dump)
            selected_frame_ref = load_system(coordx, top, frame, rsele, cpmd, dump, is_trr)
            selected_frame_ref.fix_pbc(ref_coords)
            selected_frame.fix_pbc(selected_frame_ref)
            
    if url.lower() != 'false' and url != False:
        print(f"Connecting to PyMOL at address: {url}")
        selected_frame.to_pymol(url, name)
        
    if out:
        selected_frame.write(out)
    
    print('\n \t                *****  DONE  *****                  ')
    
if __name__ == '__main__':
    main()
