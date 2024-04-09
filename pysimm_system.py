import os
import shutil
from rdkit import Chem
import pysimm
import pysimm.lmps
from pysimm import forcefield as pysimm_force_field_obj
import core_functions
import settings


class PysimmSystem:
    def __init__(
        self,
        molecule_objects_by_mol_type_nbr,
        output_dir: str = settings.output_dir,
        force_field: str = "Gaff2",
        charges: str = "Default",
    ) -> None:
        """
        Initialize a new PySIMM system instance.

        Args:
            molecule_objects_by_mol_type_nbr: Dictionary where keys are molecule types
                and items are lists of RDKit molecule objects.
            output_dir (str): Working directory for temporary files.
        """
        self.pysimm_system = None
        self.molecule_objects_by_mol_type_nbr = molecule_objects_by_mol_type_nbr
        self.output_dir = output_dir
        self._temp_mol_file = os.path.join(output_dir, "temp.mol")
        self._forcefield = getattr(pysimm_force_field_obj, force_field)()
        self._charges = charges  # Example charge method, adjust as necessary

    @staticmethod
    def _write_temp_mol_file(mol_obj, temp_mol_file: str) -> None:
        """
        Write the molecule object to a temporary MOL file.
        """
        Chem.MolToMolFile(mol_obj, temp_mol_file)
        core_functions.fix_mol_file(file=temp_mol_file)

    def _load_molecule(self, mol_file: str) -> object:
        """
        Load a molecule from a MOL file and apply a forcefield and charges.

        Args:
            mol_file (str): Path to the molecule file.

        Returns:
            The loaded molecule system.
        """
        mol_system = pysimm.system.read_mol(mol_file)
        if self._forcefield:
            mol_system.apply_forcefield(f=self._forcefield, charges=self._charges)
        return mol_system

    def initialize_system(self) -> None:
        """Create and initialize the PySIMM system with box dimensions."""
        # self.pysimm_system = pysimm.system.System()

        self._write_temp_mol_file(
            mol_obj=self.molecule_objects_by_mol_type_nbr["1"][0],
            temp_mol_file=self._temp_mol_file,
        )
        self.pysimm_system = self._load_molecule(mol_file=self._temp_mol_file)
        self.current_mol = self.pysimm_system.copy()
        self._parent_molecule = self.current_mol.copy()

        self.pysimm_system.dim.xlo = 0
        self.pysimm_system.dim.ylo = 0
        self.pysimm_system.dim.zlo = 0
        self.pysimm_system.dim.xhi = 100
        self.pysimm_system.dim.yhi = 100
        self.pysimm_system.dim.zhi = 100

    def set_box_dimensions(self, dimensions: dict) -> None:
        """Update dimensions of the PySIMM system to box_dimensions."""
        self.pysimm_system.dim.xlo = 0
        self.pysimm_system.dim.ylo = 0
        self.pysimm_system.dim.zlo = 0
        self.pysimm_system.dim.xhi = dimensions["a"]
        self.pysimm_system.dim.yhi = dimensions["b"]
        self.pysimm_system.dim.zhi = dimensions["c"]

    def load_all_molecules_into_system(self) -> None:
        """Load every molecule into the PySIMM system."""
        first_mol = True
        for mol_type, mol_list in self.molecule_objects_by_mol_type_nbr.items():
            parent_molecule = mol_list[0]  # The first molecule is the parent
            self._write_temp_mol_file(parent_molecule, self._temp_mol_file)
            parent_system = self._load_molecule(self._temp_mol_file)
            # Add the parent molecule system first
            if first_mol:
                first_mol = False
            else:
                self.pysimm_system.add(parent_system.copy(), change_dim=True)
            print(f"Currently loading {mol_type}...")

            # For every other molecule of the same type, add a copy of the parent system
            for mol in mol_list[1:]:
                self._write_temp_mol_file(mol, self._temp_mol_file)
                # Instead of reloading force field and charges, create a copy of the parent
                # and only update positions from the temp mol file.
                duplicate_system = self._load_molecule_positions_only(
                    parent_system, self._temp_mol_file
                )
                self.pysimm_system.add(duplicate_system, change_dim=True)

    @staticmethod
    def _load_molecule_positions_only(parent_system, mol_file: str) -> object:
        """
        Create a copy of the parent molecule system with updated positions from the MOL file.

        Args:
            parent_system: The parent molecule system from which to copy.
            mol_file (str): Path to the molecule file for position updates.

        Returns:
            The duplicated molecule system with updated positions.
        """
        # Create a copy of the parent system for the duplicate molecule
        duplicate_system = parent_system.copy()

        # Load the new molecule from the MOL file for position updates
        new_mol_system = pysimm.system.read_mol(mol_file)

        # Update positions in the duplicate system with those from the new molecule
        for parent_particle, new_particle in zip(
            duplicate_system.particles, new_mol_system.particles
        ):
            parent_particle.x, parent_particle.y, parent_particle.z = (
                new_particle.x,
                new_particle.y,
                new_particle.z,
            )

        return duplicate_system

    def generate_lammps_inputs(self) -> None:
        """Generate LAMMPS inputs for simulation."""
        simulation = pysimm.lmps.Simulation(self.pysimm_system, name="my_simulation")
        try:
            simulation.run(save_input=True)
        except TypeError as e:
            print(f"Warning: TypeError {e} caught and skipped.")
        shutil.move(
            src="temp.lmps", dst=os.path.join(self.output_dir, "structure.data")
        )
