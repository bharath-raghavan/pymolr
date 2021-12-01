import sys
import xmlrpc.client as xmlrpclib
import argparse
import gmx2pymol as g2p

def trr(args):
    url = 'http://localhost:9123' if args.url is None else args.url
    coords = args.trr
    name = 'mol1' if args.name is None else args.name
    top = args.top
    info = args.info
    frame = 0 if args.frame is None else int(args.frame)
    sele = 'all' if args.sele is None else ' '.join(args.sele)
    ref = args.ref
    out = args.out
    
    if top:
        a = g2p.System(mpt_file=top, trr_files=coords)
        
        if info:
            print(f"Number of frames in TRR: {a.nframes}\n")
        
        if frame == -1:
            frame = a.nframes - 1
            
        selected_frame = a.select(sele, frame=frame)
        
        if ref:
            selected_frame.fix_pbc(g2p.select_coords(top, ref, sele))
            
        if url.lower() != 'false':
            print(f"Connecting to PyMOL at address: {url}")
            selected_frame.to_pymol(url, name)
        
        if out:
            selected_frame.write(out)
    else:
        print("Topology not specified!")
        sys.exit(1)
        
def coords(args):
    url = 'http://localhost:9123' if args.url is None else args.url
    coords = args.coords
    name = 'mol1' if args.name is None else args.name
    top = args.top
    sele = 'all' if args.sele is None else ' '.join(args.sele)
    ref = args.ref

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
    
def main():
    print('\n \t                ***** gmx2pymol *****                  ')
    print('\n \t For more information type gmx2pymol [subcommand] --help \n')
    
    parser = argparse.ArgumentParser(prog='gmx2pymol')
    
    subparsers = parser.add_subparsers(title='valid subcommands', metavar='')  # Turns off list of subcommands

    #####
    parser_trr = subparsers.add_parser('trr',
                                          help='viewing and extracting frames from TRR files')
    trr_input = parser_trr.add_argument_group('options to specify input')
    
    trr_input.add_argument("trr", help="coordinate file", type=str)
    trr_input.add_argument("top", help="topology file", type=str)
    trr_input.add_argument("-i", "--info", help="print info about trr", metavar=False, type=bool)
    trr_input.add_argument("-u", "--url", help="PyMOL remote server URL", type=str, metavar='')
    trr_input.add_argument("-n", "--name", help="name of molecule in PyMOL", metavar='')
    trr_input.add_argument("-f", "--frame", help="frame to read if trr", metavar='', type=int)
    trr_input.add_argument("-r", "--ref", help="reference to fix pbc", metavar='', type=str)
    trr_input.add_argument("-s", "--sele",  nargs='+', help="selection query", metavar='')
    trr_input.add_argument("-o", "--out", help="output coordinate file", metavar='', type=str)
    trr_input.set_defaults(feature=False)
    parser_trr.set_defaults(func=trr)
    
    parser_coords = subparsers.add_parser('coords',
                                          help='viewing GRO/PDB files')
    coords_input = parser_coords.add_argument_group('options to specify input')
    
    coords_input.add_argument("coords", help="coordinate file", type=str)
    coords_input.add_argument("-t", "--top", help="topology file", type=str, metavar='')
    coords_input.add_argument("-u", "--url", help="PyMOL remote server URL", type=str, metavar='')
    coords_input.add_argument("-n", "--name", help="name of molecule in PyMOL", metavar='')
    coords_input.add_argument("-r", "--ref", help="reference to fix pbc", metavar='', type=str)
    coords_input.add_argument("-s", "--sele",  nargs='+', help="selection query", metavar='')
    coords_input.set_defaults(feature=False)
    parser_coords.set_defaults(func=coords)
    
    args = parser.parse_args()
    if vars(args) == {}:
        sys.exit()
    subcommand = args.func.__name__
    args.func(args)
      
if __name__ == '__main__':
    main()