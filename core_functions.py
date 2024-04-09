import functools
import inspect
import os
from datetime import datetime
from typing import Dict
import settings
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, Descriptors, Draw


# Global variable to track the last call signature
last_call_signature = None


def log_message_to_file(message: str) -> None:
    """
    Logs a custom message to 'calls_log.txt', prepending it with the current timestamp.

    This function retrieves the current timestamp, formats it together with the provided message,
    and appends this information to the file named 'calls_log.txt'. If the file does not exist,
    it will be created. The timestamp is formatted as 'YYYY-MM-DD_HH:MM:SS'.

    :param message: The text message to log.
    :type message: str
    """

    # Get the current timestamp
    now = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")

    # Construct the log message with timestamp
    formatted_message = f"{now} | {message}\n"

    # Open the log file and append the message
    with open("calls_log.txt", "a") as log_file:
        log_file.write(formatted_message)


def track_call_depth(func):
    """
    Decorator that logs calls to the decorated function with timestamp and call signature.
    Differentiates consecutive calls to different functions, appending details to `calls_log.txt`.
    Utilizes `inspect` for call signature and `datetime` for timestamps.

    :param func: Function to decorate.
    :return: Wrapper function with logging capability.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        global last_call_signature

        # Construct the current call signature
        frame = inspect.currentframe().f_back
        filename = os.path.basename(frame.f_code.co_filename).replace(".py", "")
        class_name = ""
        if "self" in frame.f_locals:
            class_name = frame.f_locals["self"].__class__.__name__ + "."
        current_call_signature = f"{filename}.{class_name}{func.__name__}"

        # Get the current timestamp
        now = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")

        # Check if the current call signature matches the last one
        if current_call_signature != last_call_signature:
            # Construct the log message
            log_message = f"{now} | {current_call_signature}\n"

            # Open the log file and append the log message
            with open("calls_log.txt", "a") as log_file:
                log_file.write(log_message)

            last_call_signature = current_call_signature

        # Call the original function
        return func(*args, **kwargs)

    return wrapper


def split_number(number_str: str, right_part_length: int) -> tuple:
    """
    Split a string representation of a number into left and right parts.

    :param number_str: The string representation of the number.
    :type number_str: str
    :param right_part_length: The length of the right part to extract.
    :type right_part_length: int
    :return: A tuple containing the left part and the right part of the string.
    :rtype: tuple
    """
    left_part: str = number_str[:-right_part_length]
    right_part: str = number_str[-right_part_length:]
    return left_part, right_part


@track_call_depth
def fix_mol_file(file: str) -> None:
    """
    Fix a Molecular Structure File (MOL file) by making specific modifications. The fixes focus on making sure large
    float numbers are formatted/spaced correctly (X.XX Y.XX instead of X.XXY.XX),
    which Rdkit's export function seems to have an issue with.

    :param file: The path to the MOL file to be fixed.
    :type file: str
    :return: None
    """
    nbr_of_atoms: int = 0
    nbr_of_bonds: int = 0

    # Get the full path of the input file
    full_path: str = os.path.abspath(file)
    # Get the directory name from the full path
    directory_name: str = os.path.dirname(full_path)
    # Create the full path for the output file
    write_file_path: str = os.path.join(directory_name, "fixed.mol")

    with open(file=full_path, mode="r") as read_file, open(
        write_file_path, "w"
    ) as write_file:
        for line in read_file:
            line_split = line.split()
            if len(line_split) == 10:
                left_part, right_part = split_number(
                    number_str=line_split[0], right_part_length=3
                )
                line_split[0] = " ".join([left_part, right_part])
                write_file.write(" ".join(line_split) + "\n")
            elif len(line_split) == 16:
                nbr_of_atoms += 1
                write_file.write(line)
            elif len(line_split) == 4:
                nbr_of_bonds += 1
                write_file.write(line)
            elif len(line_split) == 3:
                nbr_of_bonds += 1
                left_part, right_part = split_number(
                    number_str=line_split[0], right_part_length=3
                )
                line_split[0] = " ".join([left_part, right_part])
                write_file.write(" " + " ".join(line_split) + "\n")
            else:
                write_file.write(line)

    # Replace the original MOL file with the fixed file
    os.replace(write_file_path, full_path)


@track_call_depth
def convert_smiles_type_to_number(dictionary: dict, map_dictionary: dict) -> dict:
    """
    Converts the keys of a dictionary from SMILES to a specified numeric format.

    Iterates through the input dictionary, using a mapping dictionary to convert each key from
    a SMILES string to a corresponding numeric format, as defined in `map_dictionary`. The values
    remain unchanged.

    :param dictionary: The original dictionary with SMILES strings as keys.
    :param map_dictionary: A dictionary mapping SMILES strings to their numeric representations.
    :type dictionary: dict
    :type map_dictionary: dict
    :return: A new dictionary with keys converted to numeric format.
    :rtype: dict
    """
    new_dictionary = {}
    for key, value in dictionary.items():
        new_dictionary[map_dictionary[key]] = value
    return new_dictionary


@track_call_depth
def save_molecule_image(mol: Chem.Mol, filename: str) -> None:
    """
    Generates and saves an image representation of a molecule to a specified file.

    This function removes hydrogen atoms from the molecule, creates an image representation
    of the modified molecule, and saves it to the file specified by `filename`. The image
    is sized at 512x512 pixels and adjusted to fit the molecule representation.

    :param mol: The molecule to visualize.
    :param filename: The path and name of the file to save the image to.
    :type mol: Chem.Mol
    :type filename: str
    """

    img = Draw.MolToImage(Chem.RemoveHs(mol), size=(256, 256), fitImage=True)
    img.save(filename)


@track_call_depth
def read_pdb_file(filename: str) -> list:
    """
    Reads the contents of a PDB file and returns it as a list of lines.

    :param filename: The path to the PDB file.
    :type filename: str
    :return: A list containing the lines of the file.
    :rtype: list
    """
    with open(filename, "r") as file:
        lines = file.read().splitlines()
    return lines


@track_call_depth
def extract_molecules(pdb_lines: list) -> list:
    """
    Extracts and groups molecules from a list of PDB file lines into separate blocks.

    Parses through the given list of lines from a PDB file, identifying and grouping lines
    related to individual molecules. Each molecule block starts with the file's header and includes
    lines describing the molecule's atoms and connections. Blocks are delineated by the presence of
    "HETATM" lines for atom information and "CONECT" lines for connectivity, with a new block starting
    after each "CONECT" line. The parsing stops upon reaching the "MASTER" record, indicating the end
    of relevant molecule data.

    :param pdb_lines: The lines of a PDB file.
    :type pdb_lines: list
    :return: A list of molecule blocks, each represented as a list of lines.
    :rtype: list
    """
    pdb_molecule_blocks = []
    molecule = []
    header = pdb_lines[0]
    new_molecule = True
    for line in pdb_lines[1:]:
        if line.startswith("MASTER"):
            break
        if line.startswith("HETATM") and new_molecule:
            if molecule:
                pdb_molecule_blocks.append(molecule)
            molecule = [header]
            new_molecule = False
        molecule.append(line)
        if line.startswith("CONECT"):
            new_molecule = True
    if molecule:
        pdb_molecule_blocks.append(molecule)
    return pdb_molecule_blocks


@track_call_depth
def smiles_from_blocks(blocks: list) -> list:
    """
    Converts molecule blocks into SMILES strings with their corresponding RDKit molecule objects.

    Iterates over a list of molecule blocks, each block being a list of strings representing lines
    from a PDB file. For each block, the function constructs a molecule object using RDKit, attempts
    to determine bond orders, and generates a SMILES string representation of the molecule. The function
    returns a list of tuples, each containing a SMILES string and its corresponding RDKit molecule object.

    Note: PDB files do not contain connectivity information. Therefore, the program uses Rdkit's function
    "rdDetermineBonds.DetermineBondOrders(mol)" to deduce bonds in a molecule. The original authors of this function is
    Jan Jensen and his research group, who published "xyz2mol" to estimate bonds and bond orders from Quantum Mechanical
    simulations. Read more at https://greglandrum.github.io/rdkit-blog/posts/2022-12-18-introducing-rdDetermineBonds.html.

    :param blocks: A list of molecule blocks, each a list of strings from a PDB file.
    :type blocks: list
    :return: A list of tuples, each containing a SMILES string and an RDKit molecule object.
    :rtype: list
    """
    smiles_all = []
    for block in blocks:
        pdb_block = "\n".join(block)
        mol = Chem.MolFromPDBBlock(pdb_block, sanitize=True, removeHs=False)
        if mol:
            rdDetermineBonds.DetermineBondOrders(mol)
            smiles = Chem.MolToSmiles(mol, allBondsExplicit=True)
            smiles_all.append((smiles, mol))
    return smiles_all


@track_call_depth
def aggregate_molecules(smiles_list: list) -> tuple:
    """
    Aggregates molecules based on their SMILES representation and counts their occurrences.

    Processes a list of tuples, each containing a SMILES string and its corresponding RDKit molecule object,
    to aggregate and count occurrences of unique molecules. It returns a tuple containing two dictionaries:
    the first maps SMILES strings to a dictionary with the count of molecules and the original SMILES string;
    the second maps SMILES strings to a list of corresponding RDKit molecule objects.

    :param smiles_list: A list of tuples, each containing a SMILES string and an RDKit molecule object.
    :type smiles_list: list
    :return: A tuple containing two dictionaries:
             1. Molecule counts: Maps SMILES strings to a dictionary with the number of occurrences and the SMILES string.
             2. Molecule objects by type: Maps SMILES strings to a list of RDKit molecule objects.
    :rtype: tuple
    """
    molecule_counts = {}
    molecule_objects_by_mol_type = {}
    for n, smiles in enumerate(smiles_list):
        if smiles[0] in molecule_counts:
            molecule_counts[smiles[0]]["nbr_of_mols"] += 1
            molecule_objects_by_mol_type[smiles[0]].append(smiles[1])
        else:
            molecule_counts[smiles[0]] = {"nbr_of_mols": 1, "smiles": smiles[0]}
            # molecule_objects_by_mol_type[smiles[0]] = [smiles[1]]
            molecule_objects_by_mol_type[smiles[0]] = [smiles[1]]
    return molecule_counts, molecule_objects_by_mol_type


@track_call_depth
def smiles_to_mol_type_map(dictionary: dict) -> dict:
    """
    Maps SMILES strings to a unique numeric identifier.

    :param dictionary: A dictionary with SMILES strings as keys and dictionaries containing molecule
                       information as values.
    :type dictionary: dict
    :return: A dictionary mapping SMILES strings to unique numeric identifiers (as strings).
    :rtype: dict
    """
    n = 1
    smi_to_type_map = {}
    for key, value in dictionary.items():
        smi_to_type_map[value["smiles"]] = f"{n}"
        n += 1
    return smi_to_type_map


@track_call_depth
def compute_molecular_properties(molecule_counts: dict) -> dict:
    """
    Computes molecular properties for a set of molecules and generates their images.

    For each molecule described by its SMILES string in the input dictionary, this function calculates
    its molecular mass and chemical formula. It then generates an image of the molecule, saving it with
    a filename that includes the molecule's numeric identifier and SMILES string. The function returns a
    dictionary where each key represents a molecule identifier (e.g., "mol_1") and the value is another
    dictionary containing the molecule's chemical properties.

    :param molecule_counts: A dictionary with SMILES strings as keys and dictionaries as values,
                            where each dictionary contains the number of molecules and the SMILES string.
    :type molecule_counts: dict
    :return: A dictionary mapping molecule identifiers to their properties (SMILES string, count,
             molecular mass, and chemical formula) and generating an image for each.
    :rtype: dict
    """
    properties = {}
    for n, (smiles, info) in enumerate(molecule_counts.items(), 1):
        molecule = Chem.MolFromSmiles(info["smiles"])
        molecule_h = Chem.AddHs(molecule)
        molecular_mass = Descriptors.ExactMolWt(molecule_h)
        # Calculate the chemical formula properly
        element_counts = {}
        for atom in molecule_h.GetAtoms():
            element = atom.GetSymbol()
            element_counts[element] = element_counts.get(element, 0) + 1
        properties[f"mol_{n}"] = {
            "smiles": smiles,
            "nbr_of_mols": info["nbr_of_mols"],
            "molecular_mass": molecular_mass,
            "chemical_formula": element_counts,
        }
        # Generate and save an image for each unique molecule
        image_filename = os.path.join(
            settings.output_dir, f"mol_{n}.png"
        )  # Use a filename that includes the molecule's index and SMILES
        save_molecule_image(molecule_h, image_filename)
    return properties


@track_call_depth
def mixture_properties(properties: dict) -> dict:
    """
    Calculates the aggregate properties of a molecular mixture.

    This function computes the total mass, total number of molecules, average molecular mass,
    and combined chemical formula for a mixture based on individual molecular properties provided
    in the input dictionary.

    :param properties: A dictionary where keys are molecule identifiers and values are dictionaries
                       containing properties such as the number of molecules, molecular mass,
                       and chemical formula for each molecule type.
    :type properties: dict
    :return: A dictionary containing the total number of molecules, total mass, average molecular
             mass, and combined chemical formula of the mixture.
    :rtype: dict
    """
    total_mass = sum(
        info["nbr_of_mols"] * info["molecular_mass"] for info in properties.values()
    )
    total_molecules = sum(info["nbr_of_mols"] for info in properties.values())
    # Initialize an empty dictionary for the chemical formula of the mixture
    chemical_formula = {}
    for info in properties.values():
        for element, count in info["chemical_formula"].items():
            if element in chemical_formula:
                chemical_formula[element] += count * info["nbr_of_mols"]
            else:
                chemical_formula[element] = count * info["nbr_of_mols"]
    return {
        "nbr_of_molecules": total_molecules,
        "total_mass": total_mass,
        "avg_molecular_mass": total_mass / total_molecules,
        "chemical_formula": chemical_formula,
    }


@track_call_depth
def extract_pdb_box_dimensions(pdb_header: str) -> Dict[str, float]:
    """
    Extracts box dimensions from PDB header.

    Parameters:
    - pdb_header (str): PDB header containing box dimensions.

    Returns:
    - dimensions (Dict[str, float]): Dictionary containing box dimensions.
    """
    split_line = pdb_header.split()
    return {
        "a": float(split_line[1]),
        "b": float(split_line[2]),
        "c": float(split_line[3]),
        "box.alpha": float(split_line[4]),
        "box.beta": float(split_line[5]),
        "box.gamma": float(split_line[6]),
    }
