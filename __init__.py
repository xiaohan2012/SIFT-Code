__author__ = 'xiaohan'
from schrodinger import structure
import fp_gen1
import avg_sift
import schrodinger.utils.fileutils as fileutils
import os


def start_with_atom(line):
    """
    tell if the line is atom line
    """
    return line.split()[0].lower() == 'atom'

def belongs_to_antibody(line):
    """
    tell if the line belongs to antibody
    """
    try:
        return line.split()[4].upper() in 'ABHL'
    except:
        return False


def split_protein_protein_complex_manual(st_path,root_dir):
    """
    @params:complex structure file path
    @return antibody structure, antigen structure
    """
    complex_id,_ = fileutils.splitext(os.path.basename(st_path))

    if not os.path.exists(root_dir):
        os.mkdir(root_dir)
    result_save_path = os.path.join(root_dir,complex_id)
    if not os.path.exists(result_save_path):
        os.mkdir(result_save_path)

    ab_fname = '%s/%s_antibody.pdb' %(result_save_path,complex_id)
    ag_fname = '%s/%s_antigen.pdb' %(result_save_path,complex_id)
    with open(ab_fname,'w') as ab_f:
        with open( ag_fname,'w') as ag_f:
            with open(st_path) as f:
                for line in f.readlines():
                    if start_with_atom(line):
                        if belongs_to_antibody(line):
                            ab_f.write(line)
                        else:
                            ag_f.write(line)

    antibody_id,_ = fileutils.splitext(os.path.basename(ab_fname))
    antigen_id,_ = fileutils.splitext(os.path.basename(ag_fname))

    antibody = structure.StructureReader(ab_fname).next()
    antibody.title = antibody_id

    antigen = structure.StructureReader(ag_fname).next()
    antigen.title = antigen_id

    return antibody,antigen

def belongs_to_receptor(line):
    aa = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    try:
        return line.split()[3] in aa
    except:
        return False

def belongs_to_molecule(line):
    try:
        return line.split()[0] == 'HETATM'
    except :
        return False

def split_protein_molecule_complex_manual(st_path,root_dir = 'processed_data_protein_molecule'):
    """
    @params:complex structure file path
    @return receptor,small molecule
    """
    complex_id,_ = fileutils.splitext(os.path.basename(st_path))
    root_dir
    if not os.path.exists(root_dir):
        os.mkdir(root_dir)

    result_save_path = os.path.join(root_dir,complex_id)

    if not os.path.exists(result_save_path):
        os.mkdir(result_save_path)

    rec_fname = '%s/%s_receptor.pdb' %(result_save_path,complex_id)
    bind_fname = '%s/%s_binder.pdb' %(result_save_path,complex_id)
    with open(rec_fname,'w') as rec:
        with open( bind_fname,'w') as bind:
            with open(st_path) as f:
                for line in f.readlines():
                    if start_with_atom(line):
                        rec.write(line)
                    elif belongs_to_molecule(line):
                        bind.write(line)
    rec_id,_ = fileutils.splitext(os.path.basename(rec_fname))
    bind_id,_ = fileutils.splitext(os.path.basename(bind_fname))

    rec = structure.StructureReader(rec_fname).next()
    rec.title = rec_id

    bind = structure.StructureReader(bind_fname).next()
    bind.title = bind_id

    return rec,bind
    
def split_protein_protein_complex(st_path):
    st = structure.StructureReader(st_path).next()
    chains = [c for c in st.chain]
    print chains

    antibody = structure.create_new_structure()

    
    antibody.title =  '%s_antibody' %(st.title)
    antigen = structure.create_new_structure()
    antigen.title = '%s_antigen' %(st.title)
    for chain in chains:
        # considering this chain

        for a in st.chain[chain.name].atom:
            #if a.bond_total == 0:continue
            if chain.name in 'HLAB':# it is antibody
                antibody.addAtom(a.element,a.x,a.y,a.z,a.color,a.atom_type)
            else:
                antigen.addAtom(a.element,a.x,a.y,a.z,a.color,a.atom_type)

    antibody.write('%s.pdb' %antibody.title)
    antigen.write('%s.pdb' %antigen.title)
    return antigen, antibody

def gen_protein_protein_complex_avg_sift(complex_st_path,processed_data_path = 'processed_data_protein_protein'):
    complex_id,  ext = fileutils.splitext(os.path.basename(complex_st_path))

    sift_path = '%s/%s/%s_pattern.dat' %(processed_data_path,complex_id,complex_id)

    if os.path.exists(sift_path):
        print complex_id,'is processed'
        return


    antibody, antigen = split_protein_protein_complex_manual(complex_st_path,processed_data_path)#load complex structure and split it

    fp_gen_path = fp_gen1.gen_fp(receptor=antibody,binder = antigen,complex_id = complex_id,root_path = processed_data_path)#get the finger print

    avg_sift.gen_avg_sift(fp_gen_path,sift_path)#generate average sift

def gen_protein_molecule_complex_avg_sift(complex_st_path):
    complex_id,  ext = fileutils.splitext(os.path.basename(complex_st_path))


    sift_path = '%s/%s/%s_pattern.dat' %(processed_data_path,complex_id,complex_id)
    
    if os.path.exists(sift_path):
        print complex_id,'is processed'
        return
    
    rec, bind = split_protein_molecule_complex_manual(complex_st_path)#load complex structure and split it

    fp_gen_path = fp_gen1.gen_fp(receptor=rec,binder = bind,complex_id = complex_id,root_path = processed_data_path)#get the finger print

    avg_sift.gen_avg_sift(fp_gen_path,sift_path)#generate average sift

def process_test_data(protein_path = 'test_protein',ligand_path = 'test_ligand',o_path = 'test_output'):
    for ligand_file in os.listdir(ligand_path):
        
        ligand_st_path = os.path.join(ligand_path,ligand_file)#Ligand file path
        protein_st_path = os.path.join(protein_path,ligand_file)#protein file path
        
        print ligand_st_path,protein_st_path


        if not os.path.exists(protein_st_path):#protein file not exist
            #print protein_st_path,'not exist'
            continue
        else:
            ligand_st = structure.StructureReader(ligand_st_path).next()
            protein_st = structure.StructureReader(protein_st_path).next()

            complex_id=ligand_file[:4]
            try:#dangerous
                fp_gen_path = fp_gen1.gen_fp(receptor=protein_st,binder=ligand_st,complex_id=complex_id,root_path=o_path)
            except:
                print complex_id,'cannot be processed'
            sift_path = os.path.join(complex_id,'pattern.dat')
            if os.path.exists(sift_path):
                print complex_id,'is processed'
                continue
            avg_sift_path = os.path.join(o_path,sift_path)
            avg_sift.gen_avg_sift(fp_gen_path,avg_sift_path)#generate average sift

def get_avg_matrix(receptor_path='',ligands_path = [],output_fp_path = '',output_pattern_path = ''):
    #get receptor
    receptor = structure.StructureReader(receptor_path).next()
    receptor.title = receptor_path
    print receptor_path,'set'

    rec_tree = fp_gen1.distance_tree()

    rec_tree.set_receptor_structure(receptor)
    rec_tree.parse_receptor()

    #get ligands
    ligands = []
    for ligand_path in ligands_path:
        lig = structure.StructureReader(ligand_path).next()
        lig.title = ligand_path
        ligands.append(lig)

        print ligand_path,'added'

    for lig in ligands:
        rec_tree.find_close_residues(lig)
        print lig,'close residues found'


    with open(output_fp_path,  'w') as out_fp:
        for key in rec_tree.fingerprints.sifts.keys():
            rec_tree.fingerprints.fill_missing_zeros(rec_tree.min_res,  rec_tree.max_res,  key)
            fp_string = rec_tree.fingerprints.get_sift_string(key)
            out_fp.write(rec_tree.receptor.title + ':' + key + ':' + str(rec_tree.min_res) + ':' + fp_string + '\n')
        print 'finger print',output_fp_path,'saved'

    avg_sift.gen_avg_sift(output_fp_path,output_pattern_path)
    
    print 'pattern',output_pattern_path,'saved'

if __name__ == '__main__':

    """
    data_path = 'protein2molecule'
    
    for pdb_file in os.listdir(data_path):
        try:
            gen_protein_molecule_complex_avg_sift(os.path.join(data_path,pdb_file))
            #split_protein_molecule_complex_manual(os.path.join(data_path,pdb_file))
            print pdb_file,'ok'
        except:
            print pdb_file,'failed'
    """

    """
    data_path = 'protein-protein-data'
    for pdb_file in os.listdir(data_path):
        #if pdb_file in ['1BJ1.pdb','1JRH.pdb','1CE1.pdb','1JTO.pdb','3HFM']:continue
        try:
            complex_st_path = '%s/%s' %(data_path,pdb_file)
        #split_protein_protein_complex_manual(complex_st_path)
            gen_protein_protein_complex_avg_sift(complex_st_path)
        except :
            print pdb_file,'failed'

    """
    """
    data_dir = 'processed_data_protein_protein/1HEZ'
    receptor_file = 'abcd'
    receptor_path = os.path.join(data_dir,receptor_file +'.pdb')

    ligand_file = 'e'
    ligand_files = [ligand_file +'.pdb']

    ligands_path = map(lambda s:os.path.join(data_dir,s),ligand_files)

    output_fp_path = os.path.join(data_dir,'%s_%s_fp.dat' %(receptor_file,ligand_file))
    output_pattern_path = os.path.join(data_dir,'%s_%s_pattern.dat' %(receptor_file,ligand_file))
    
    get_avg_matrix(receptor_path,ligands_path,output_fp_path,output_pattern_path)
"""
    gen_protein_protein_complex_avg_sift('1918_complex.pdb')
