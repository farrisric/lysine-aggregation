from ase import Atoms
from ase.io import write,read
import re 

class ChemicalShift():
    def __init__(self, nmr_out):
        self.nmr_out = None
        self.atoms_obj = None
        self.symbols = []
        self.positions = []
        self.chemical_shifts = []
        
        self._open_nmr_out(nmr_out)

    
    def _open_nmr_out(self, nmr_out):
        try:
            self.atoms_obj = read(nmr_out)
            self.symbols = self.atoms_obj.symbols
            self.positions = self.atoms_obj.positions
            with open(nmr_out) as nmr_out:
                nmr_out = nmr_out.readlines()
                self.nmr_out = nmr_out
        except:
            with open(nmr_out) as nmr_out:
                nmr_out = nmr_out.readlines()
                self.nmr_out = nmr_out
      

    def _get_positions_and_symbols_from_nmr_out(self):
        reg = ' [OCNH],(([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[eE]([+-]?\d+))?(,([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[eE]([+-]?\d+))?)+)'
        
        for line in self.nmr_out:
            if not len(line) > 0:
                continue
            if re.search(reg, line) and len(line.split(',')) == 5:  
                symbol, _, x, y ,z = line.split(',')
                self.symbols.append(symbol.replace(' ', ''))
                self.positions.append((x, y ,z))

    def get_chemical_shift_from_nmr_out(self):
        prefactors = {'C' : 161.5794+24.637, 'N': 213.8677+70, 'H' : 31.78218333+0, 'O' : 0}

        reg = "Anisotropy"
        for line in self.nmr_out:
            if not len(line) > 0:
                continue
            if re.search(reg, line):
                _, symbol, _,_, anisotropy, _,_,_ = line.split()
                if symbol == 'O' or symbol == 'Si':
                    self.chemical_shifts.append(0)
                else:
                    self.chemical_shifts.append(prefactors[symbol] - float(anisotropy))

    def get_atoms_object(self):
        if len(self.symbols) == len(self.chemical_shifts):
            return Atoms(
                symbols=self.symbols,
                positions=self.positions,
                charges=self.chemical_shifts)
        else:
            return False