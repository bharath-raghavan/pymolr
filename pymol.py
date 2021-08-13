import sys
import xmlrpc.client as xmlrpclib
import argparse
import gmx2pymol as g2p

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("coords", help="coordinate file", type=str)
    parser.add_argument("-t", "--top", help="topology file", type=str, metavar='')
    parser.add_argument("-u", "--url", help="PyMOL remote server URL", type=str, metavar='')
    parser.add_argument("-n", "--name", help="name of molecule in PyMOL", metavar='')
    parser.add_argument("-f", "--frame", help="frame to read if trr", metavar='', type=int)
    parser.add_argument("-r", "--ref", help="reference to fix pbc", metavar='', type=str)
    parser.add_argument("-s", "--sele",  nargs='+', help="selection query", metavar='')
    parser.add_argument("-o", "--out", help="output coordinate file", metavar='', type=str)
    parser.set_defaults(feature=False)
    args = parser.parse_args()
    
    url = 'http://localhost:9123' if args.url is None else args.url
    coords = args.coords
    name = 'mol1' if args.name is None else args.name
    top = args.top
    frame = 0 if args.frame is None else int(args.frame)
    sele = 'all' if args.sele is None else ' '.join(args.sele)
    ref = args.ref
    out = args.out
    
    if coords.split('.')[-1] != 'trr':
        if sele == 'all':
            pymol = xmlrpclib.ServerProxy(url)
            with open(coords, 'r') as f:
                s = f.read().replace('\n','@')
            pymol.do(f'g2p_load("{s}", "{name}")')
        else:
            if not top:
                print("Topology not specified!")
                sys.exit(1)
            else:
                g2p.select_coords(top, coords, sele).to_pymol(url , name)
    else:
        if top:
            a = g2p.System(mpt_file=top, trr_files=coords)
            print(a.nframes)
            selected_frame = a.select(sele, frame=frame)
            if ref:
                selected_frame.fix_pbc(g2p.select_coords(top, ref, sele))
            selected_frame.to_pymol(url, name)
            
            if out:
                
                selected_frame.write(out)
        else:
            print("Topology not specified!")
            sys.exit(1)
        
        
if __name__ == '__main__':
    main()