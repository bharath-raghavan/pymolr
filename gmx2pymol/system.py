from os import path
import xmlrpc.client as xmlrpclib
import pickle
import numpy as np
import pandas as pd
from mimicpy import Mpt, CoordsIO
from .trr import read_trr, get_trr_frames, write_trr, TRRReader
from .rmsd import kabsch_rmsd, kabsch_weighted_fit

def select_coords(mpt, coords, selection):
    sele = Mpt.from_file(mpt).select(selection)
    ids = sele.index - 1
    coords_handle = CoordsIO(coords, 'r')
    box = np.diag(coords_handle.box)

    x = coords_handle.coords[['x','y','z']].to_numpy()[ids]
    try:
        v = coords_handle.coords[['v_x','v_y','v_z']].to_numpy()[ids]
    except KeyError:
        v = np.array([])

    return SelectedFrame(0, 0, box, sele, x, v, np.array([]), None)


class System:
    def __init__(self, mpt_file, trr_files, client=None, status_file=None, status=None):
        self.mpt = Mpt.from_file(mpt_file)
        
        if not isinstance(trr_files, list):
            trr_files = [trr_files]
            
        each_trr_frames = [get_trr_frames(trr) for trr in trr_files]
        self._trrs = list(zip(trr_files, each_trr_frames))
        
        self.nframes = sum(each_trr_frames)
        
        self.client = client
        
    @staticmethod
    def __get_frame(sele, trr_files, frame):
        if isinstance(sele, list) and len(sele) > 1:
            return [System.__get_frame(s, trr_files, frame) for s in sele]
        elif isinstance(sele, list):
            sele = sele[0]
        
        ids = sele.index - 1

        frames_so_far = 0
        trr_file_names, nframes = list(zip(*trr_files))
        trr_data = None
        for i, f in enumerate(nframes):
            if frame < frames_so_far+f:
                header, data = read_trr(trr_file_names[i], frame-frames_so_far)
                trr_data = {**header, **data}
                break
            else:
                frames_so_far += f
        
        if trr_data is None:
            raise Exception("Requested frame is greater than total number of frames.")
        
        x = trr_data['x'][ids] if 'x' in trr_data else np.array([])
        v = trr_data['v'][ids] if 'v' in trr_data else np.array([])
        f = trr_data['f'][ids] if 'f' in trr_data else np.array([])
        
        box = [trr_data['box'][0,0], trr_data['box'][1,1], trr_data['box'][2,2]] # only rect box supported for now
        
        return SelectedFrame(trr_data['step'], trr_data['time'], box, sele, x, v, f, header)
    
    @staticmethod
    def log_pickle(x, file_name):
        if file_name is None: return x
        
        xyz = np.array(x)
        
        try:
            if path.isfile(file_name):
                with open(file_name, 'rb') as f:
                    try:
                        old = pickle.load(f)
                        xyz = np.vstack((old, xyz))
                    except EOFError:
                        pass

            with open(file_name, 'wb') as f: pickle.dump(xyz, f)
        except NameError:
            pass
        
        return x
        
    def select(self, *selection, calc=None, frames=[], frame=0, as_futures=False, **kwargs):
        if isinstance(selection, str):
            sele = self.mpt.select(selection)
        else:
            sele = [self.mpt.select(s) for s in selection]
        
        trr_file = self._trrs
        
        if calc:
            do = lambda frame, **kwargs: calc(System.__get_frame(sele, trr_file, frame), **kwargs)
        else:
            do = lambda frame: System.__get_frame(sele, trr_file, frame)
        
        if frames == []:
            ret = do(frame)
        elif self.client is None:
            ret = [do(i) for i in frames]
        elif as_futures:
            ret = [self.client.submit(do, i, **kwargs) for i in frames]
        else:
            ret = [self.client.submit(do, i, **kwargs).result() for i in frames]
        
        return ret
    
    def select_coords(self, selection, filename):
        sele = self.mpt.select(selection)
        ids = sele.index - 1
        coords_handle = CoordsIO(filename, 'r')
        box = np.diag(coords_handle.box)
        
        x = coords_handle.coords[['x','y','z']].to_numpy()[ids]
        try:
            v = coords_handle.coords[['v_x','v_y','v_z']].to_numpy()[ids]
        except KeyError:
            v = np.array([])
            
        return SelectedFrame(0, 0, box, sele, x, v, np.array([]), None)
    

class SelectedFrame:
    
    def __init__(self, step, time, box, df, x, v, f, header):
        self.__header = header
        self.step = step
        self.time = time
        self.box = box
        
        self.natoms = len(df)
        if self.natoms == 1:
            one_only = True
        else:
            one_only = False
        
        if one_only:
            self.ids = df.index.to_numpy()[0]
        else:
            self.ids = df.index.to_numpy()
        
        self.__repr_list = []
        
        for column in df:
            if column == 'resid': dtype = np.int32
            elif column in ['charge', 'mass']: dtype = np.float32
            else: dtype = np.str
            self.__repr_list.append(column)
            
            if one_only:
                setattr(self, column, df[column].to_numpy(dtype=dtype)[0])
            else:
                setattr(self, column, df[column].to_numpy(dtype=dtype))
        
        if one_only and x != []:
            self.positions = x[0]
        else:
            self.positions = x
        
        if one_only and v != []:
            self.velocities = v[0]
        else:
            self.velocities = v
        
        if one_only and f != []:
            self.forces = f[0]
        else:
            self.forces = f
        self.i = 0
        self.__pbc_box = None
    
    def __getitem__(self, k):
        if isinstance(k, int): k = [k]
        dct = {}
        for i in self.__repr_list:
            dct[i] = getattr(self, i)[k]
        mpt_df = pd.DataFrame(dct, index=self.ids[k])
        
        try:
            pos = self.positions[k]
        except (TypeError,IndexError):
            pos = []
        
        try:
            v = self.velocities[k]
        except (TypeError,IndexError):
            v = []
        
        try:
            f = self.forces[k]
        except (TypeError,IndexError):
            f = []
            
        return SelectedFrame(self.step, self.time, self.box, mpt_df, pos, v, f, self.__header)
    
    def sort(self, by=None, in_place=True):
        if by is None:
            by = self.ids
        sorted_ids = np.argsort(by)
        frame = self.__getitem__(sorted_ids)
        if in_place:
            self.__dict__.update(frame.__dict__)
        else:
            return frame
        
    def __repr__(self):
        s = f"\tStep: {self.step}\tTime: {self.time}\t#Atoms: {self.natoms}\n==========================================================\n\n"

        for k in self._SelectedFrame__repr_list:
            v = getattr(self, k)
            s += f"{k}: {v}\n"

        if len(self.positions) != 0:
            s += f"\nPositions:\n==========\n{self.positions}\n"

        if len(self.velocities) != 0:
            s += f"\nVelocities:\n==========\n{self.velocities}\n"

        if len(self.forces) != 0:
            s += f"\nForces:\n==========\n{self.forces}\n"
        
        return s
    
    def __add__(self, other):
        dct = {}
        for i in self.__repr_list:
            dct[i] = np.append(getattr(self, i), getattr(other, i))
        mpt_df = pd.DataFrame(dct, index=np.append(self.ids, other.ids))

        pos = np.vstack((self.positions, other.positions))
        v = np.vstack((self.velocities, other.velocities))
        f = np.vstack((self.forces, other.forces))
                              
        return SelectedFrame(self.step, self.time, self.box, mpt_df, pos, v, f, self.__header)
    
    def __iter__(self):
        return self

    def __next__(self):
        num = self.i
        self.i += 1
        if self.i > self.natoms:
            self.i = 0
            raise StopIteration
            
        return self.__getitem__(num)
    
    def get_around(self, sele, dist=3, inplace=True, fix_pbc=True):
        if fix_pbc and self.__pbc_box:
            sele.fix_pbc(self)
        
        if sele.positions.shape != (3,):
            cen = np.mean(self.positions, axis=0)
            
        ids = np.where( np.linalg.norm(sele.positions - cen, axis=1) <= dist)[0].tolist()
        if inplace:
            return (self + sele[ids])
        else:
            return sele[ids]

    def where(self, bool_id):
        return self.__getitem__(np.where(bool_id)[0].tolist())
    
    def write(self, file, extension='gro', append=False):
        if file:
            extension = file.split('.')[-1]
        if extension == 'trr' and self.__header:
            data = self.__header
            data['box'] = np.array([[self.box[0], 0, 0], [0, self.box[1], 0], [0, 0, self.box[1]]])
            data['x'] = self.positions
            #if len(self.velocities) != 0:
            #    data['v'] = self.velocities
            #if len(self.forces) != 0:
            #    data['f'] = self.forces
            write_trr(file, data, data['endian'], data['double'], append)
        else:
            dct = {}
            for i in self.__repr_list:
                dct[i] = getattr(self, i)
            mpt_df = pd.DataFrame(dct)
            mpt_df['id'] = self.ids

            dct = {}
            dct['x'] = self.positions[:, 0]
            dct['y'] = self.positions[:, 1]
            dct['z'] = self.positions[:, 2]
            coords_df = pd.DataFrame(dct)
            coords_df['id'] = self.ids

            if file is None:
                return CoordsIO('dummy.'+extension, mode='w').write(mpt_df, coords_df, as_str=True)
            else:
                CoordsIO(file, mode='w').write(mpt_df, coords_df)

    def to_pymol(self, url='http://localhost:9123', name='mol1'):
        pymol = xmlrpclib.ServerProxy(url)
        pymol.zoom()
        
        s = self.write(None).replace('\n','@')
        
        try:
            pymol.do(f'g2p_load("{s}", "{name}")')
        except Fault:
            raise Exception("The pymolrc.py file included with gmx2pymol is not sourced on the client side")
            
    def correct_box(self):
        if self.positions != []:
            self.positions -= np.min(self.positions, axis=0) # make origin zero
            self.box = np.max(self.positions, axis=0)
        else:
            raise Exception
    
    def rmsd(self, other):
        pos1 = self.positions
        pos2 = other.positions
        if pos1.shape != pos2.shape:
            raise Exception

        return kabsch_rmsd(pos1, pos2, self.mass, True)
    
    def fit(self, other):
        pos1 = self.positions
        pos2 = other.positions
        if pos1.shape != pos2.shape:
            raise Exception
        
        return kabsch_weighted_fit(pos1, pos2, self.mass, True)
    
    def __get_subset(self, s, frame, dims):
            box = [(dims[0],0), (dims[1],0), (dims[2],0)]

            tup = (int(s[0]), int(s[1]), int(s[2]))

            for i, t in enumerate(tup):
                if t == 1:
                    box[i] = (dims[i]/2, -dims[i]/2)
            
            return frame.where(
                (frame.positions[:, 0] < box[0][0]) & 
                (frame.positions[:, 0] >= box[0][1]) & 
                (frame.positions[:, 1] < box[1][0]) & 
                (frame.positions[:, 1] >= box[1][1]) & 
                (frame.positions[:, 2] < box[2][0]) & 
                (frame.positions[:, 2] >= box[2][1])).sort(in_place=False)
        
    def __get_images(self, box):
        b1, b2, b3 = box[0], box[1], box[2]
        
        import copy
        sele1 = copy.deepcopy(self)
        sele2 = copy.deepcopy(self)
        sele3 = copy.deepcopy(self)
        sele4 = copy.deepcopy(self)
        sele5 = copy.deepcopy(self)
        sele6 = copy.deepcopy(self)
        sele7 = copy.deepcopy(self)
        
        sele1.positions[:,0] -= b1
        sele2.positions[:,1] -= b2
        sele3.positions[:,2] -= b3

        sele4.positions[:,0] -= b1
        sele4.positions[:,1] -= b2

        sele5.positions[:,0] -= b1
        sele5.positions[:,2] -= b3

        sele6.positions[:,0] -= b1
        sele6.positions[:,1] -= b2
        sele6.positions[:,2] -= b3

        sele7.positions[:,1] -= b2
        sele7.positions[:,2] -= b3
        
        return self+sele1+sele2+sele3+sele4+sele5+sele6+sele7
   
        
    def fix_pbc(self, ref, in_place=True):
        
        self.positions -= np.min(self.positions, axis=0) # make origin zero
        box = np.max(self.positions, axis=0)
            
        sel = self.__get_images(box)
        
        if isinstance(ref, str):
            pbc_box = ref
        else:
            pbc_box = ref.__pbc_box
        
        if pbc_box is None:
            if ref.natoms != self.natoms:
                raise Exception("Ref not same as self")
            
            rmsd = {}
            for i in range(8):
                s = "{0:03b}".format(i)
                e = self.__get_subset(s, sel, box)
                if e.natoms == self.natoms:
                    pos, rmsd_val = e.fit(ref)
                    e.positions = pos
                    rmsd[s] = [e, rmsd_val]
                else:
                    raise Exception(f"{s} {e.natoms}")
            
            min_rmsd = min(rmsd, key=lambda x: rmsd[x][1])
            frame = rmsd[min_rmsd][0]
            frame.__pbc_box = min_rmsd
        else:
            frame = self.__get_subset(pbc_box, sel, box)
            frame.__pbc_box = pbc_box
        
        if in_place:
            self.__dict__.update(frame.__dict__)
        else:
            return frame