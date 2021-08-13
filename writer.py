import sys
import xmlrpc.client as xmlrpclib
import argparse
import gmx2pymol as g2p

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("coords", help="coordinate file", type=str)
    parser.add_argument("-t", "--top", help="topology file", type=str, metavar='')
    parser.add_argument("-f", "--frame", help="frame to read if trr", metavar='', type=int)
    parser.add_argument("-s", "--sele",  nargs='+', help="selection query", metavar='')
    parser.set_defaults(feature=False)
    args = parser.parse_args()
    
    coords = args.coords
    top = args.top
    frame = 0 if args.frame is None else int(args.frame)
    sele = 'all' if args.sele is None else ' '.join(args.sele)
    
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
                wfg.select_coords(top, coords, sele).to_pymol(url , name)
    else:
        if top:
            a = g2p.System(mpt_file=top, trr_files=coords)
            a.select(sele, frame=frame).to_pymol(url, name)
        else:
            print("Topology not specified!")
            sys.exit(1)
        
        
if __name__ == '__main__':
    main()