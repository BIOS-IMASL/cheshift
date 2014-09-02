'''
cheshift plugin-in: Validate your protein model with PyMOL
Described at PyMOL wiki: http://www.pymolwiki.org/index.php/cheshift

Author : Osvaldo Martin
email: aloctavodia@gmail.com
Date: June 2014
License: GNU General Public License
Version 3.5
'''

import Tkinter
from Tkinter import *
import tkFileDialog
import Pmw
from pymol import cmd, stored
import sys
import re
import os

try:
    import numpy as np
except ImportError:
    print '<'*80 + '\n\nCheShift needs NumPy to be installed in your system, please follow the instructions at\nhttp://www.pymolwiki.org/index.php/cheshift\n\n' + '>'*80

try:
    from scipy.interpolate import griddata
except ImportError:
    print '<'*80 + '\n\nCheShift needs SciPy to be installed in your system, please follow the instructions at\nhttp://www.pymolwiki.org/index.php/cheshift\n\n' + '>'*80


path = os.path.dirname(__file__)

def __init__(self):
    """Add this Plugin to the PyMOL menu"""
    self.menuBar.addmenuitem('Plugin', 'command',
                            'Cheshift',
                            label = 'Cheshift',
                            command = lambda : mainDialog())


def validation(pdb_filename, cs_path):
    """Run the CheShift validation routine"""
    cmd.set('suspend_updates', 'on')
    cs_exp = os.path.split(cs_path)[1]
    cs_exp_name = cs_exp.split('.')[0]
    pose, residues, total_residues, states = pose_from_pdb(pdb_filename)
    reference = bmrb2cheshift(cs_exp, cs_exp_name)
    ok, ocslist_full_new = check_seq(residues, cs_exp_name)
    if  ok == 1:
        Db = load(path)
        cs_2_colors(cs_exp_name, pose, residues, total_residues, states, reference, Db)
        clean(pose)
        colorize()
    cmd.set('suspend_updates', 'off')
        #print '<'*80 + '\nCheShift-2 Validation Report saved at\n%s.zip\n' % details_path + '>'*80


def prediction(pdb_filename):
    """Run the CheShift CS prediction routine"""
    cmd.set('suspend_updates', 'on')
    pose, residues, total_residues, states = pose_from_pdb(pdb_filename)
    Db = load(path)
    raw(pose, residues, total_residues, states, Db) 
    print '<'*80 + '\nYou didn`t provide a file with chemical Shifts, hence CheShift-2 assumed you\n only wanted the predicted CS. The predicted chemical shifts can be found in the file %s.txt\n' % pose + '>'*80
    for sel in ['A', 'B', 'C', 'D']:
        cmd.delete(sel)
    cmd.set('suspend_updates', 'off')


def mainDialog():
    """ Creates the GUI """
    master = Tk()
    master.title(' CheShift ')
    w = Tkinter.Label(master, text='\nCheShift: Validate your protein model with PyMOL',
                                background = 'black',
                                foreground = 'white')
    w.pack(expand=1, fill = 'both', padx=4, pady=4)
############################ NoteBook #########################################
    Pmw.initialise()
    nb = Pmw.NoteBook(master, hull_width=430, hull_height=250)
    p1 = nb.add('Run CheShift')
    p2 = nb.add(' Color code ')
    p4 = nb.add('    About   ')
    nb.pack(padx=5, pady=5, fill=BOTH, expand=1)
############################ RUN CheShift TAB #################################
# select files
    group = Pmw.Group(p1, tag_text='Select your file')
    group.pack(fill='both', expand=1, padx=5, pady=5)
    Label(group.interior(), text =u"""
If you don't select a file, CheShift will predict the 
13C\u03B1 and 13C\u03B2 chemical shifts values for the 
currently loaded protein.
If you choose a file, CheShift will validate your protein model.
""",justify=LEFT).pack()
    Button(group.interior(), text='Chemical Shift file', command=retrieve_cs).pack()
# Run
    Button(p1, text="Run", command=run).pack(side=BOTTOM)
############################ COLOR TAB ########################################
    Label(p2, text =u"""
Colors indicate the difference between predicted and
observed 13C\u03B1 and 13C\u03B2 chemical shifts values
averaged over all uploaded conformers.
Green, yellow and red colors represent small, medium
and large differences, respectively. 
White is used if either the prediction fail or the
observed value is missing
CheShift-2 provied alternative rotamers for blue residues.
""",justify=LEFT).pack()
    Button(p2, text="Reset View", command=colorize).pack(side=BOTTOM)
############################ About TAB ########################################
    Label(p4, text = """
If you find CheShift useful please cite:

Martin O.A. Arnautova Y.A. Icazatti A.A. 
Scheraga H.A. and Vila J.A. 
A Physics-Based Method to Validate and 
Repair Flaws in Protein Structures. 
Proc Natl Acad Sci USA 2013. 110(42):16826-31
""",justify=CENTER).pack()
    master.mainloop()


def retrieve_cs():
    """Loads a Chemical Shift file provided by the user"""
    global cs_path
    cs_path = tkFileDialog.askopenfilename(title="Open chemical shift file", filetypes=[("All files","*")])
    if len(cs_path) == 0:
        del cs_path


def colorize():
    """Color the protein according to the b-factor. Uses the 
    'cheshift-color-code''"""
    try:
        cmd.spectrum('b', 'red_yellow_green', minimum='-1.0', maximum='0.0')
        cmd.select('missing', 'b = -2.0')
        cmd.color('white','missing')
        cmd.delete('missing')
        cmd.select('fixable', 'b = 2.0')
        cmd.color('blue','fixable')
        cmd.delete('fixable')
        cmd.hide()
        cmd.show('cartoon')
    except:
        pass


def run():
    """Checks if files were provided and calls the validation 
    or prediction routine"""
    pdb = 0
    cs = 0
    try:
        pdb_filename = cmd.get_names('all')[0]
        print pdb_filename
        pdb = 1
    except:
        Pmw.MessageDialog(title = 'Error',message_text = 'Please choose a\n PDB file')
    if pdb == 1:
        try:
            cs_path
            cs = 1
        except:
                pass
    if pdb == 1 and cs == 1:
        validation(pdb_filename, cs_path)
    if pdb == 1 and cs == 0:
        prediction(pdb_filename)


def bmrb2cheshift(cs_exp, cs_exp_name):
    """Parse the experimental chemical shifts file. Stores the data in an easy 
    format for further processing"""
    for line in open('%s' % cs_exp).readlines():
        if 'DSS' in line:
            reference = 1.7
        elif 'TSP' in line:
            reference = 1.82
        elif 'TMS' in line:
            reference = 0.00
    try:
        reference
    except:
        print 'Unknow reference value, please check the BMRB file'
        sys.exit()
    try:
        cs_exp_ca = []
        cs_exp_cb = []
        a = re.compile('[0-9]{1,5}\s{1,4}[A-Z]{3}\sCA\s{0,3}C.{0,5}[0-9]*\.[0-9]{0,2}')
        b = re.compile('[0-9]{1,5}\s{1,4}[A-Z]{3}\sCB\s{0,3}C.{0,5}[0-9]*\.[0-9]{0,2}')
        for line in open('%s' % cs_exp).readlines():
            if a.search(line):
                data = a.search(line).group().split()
                cs_exp_ca.append(data)
            if b.search(line):
                data = b.search(line).group().split()
                cs_exp_cb.append(data)
        len_a = len(cs_exp_ca)
        len_b = len(cs_exp_cb)
        if len_a > len_b:
            dif = len_a - len_b
            for i in range(0, dif):
                cs_exp_cb.append(['99999'])
        elif len_a < len_b:
            dif = len_b - len_a
            for i in range(0, dif):
                cs_exp_ca.append(['99999'])
        count_ca = 0
        count_cb = 0
        ocs_list = []
        while True:
            try:
                resn_ca = int(cs_exp_ca[count_ca][0])
                resn_cb = int(cs_exp_cb[count_cb][0])
                if resn_ca == resn_cb:
                    line = '%4s %3s  %6.2f  %6.2f\n' % (cs_exp_ca[count_ca][0], cs_exp_ca[count_ca][1], float(cs_exp_ca[count_ca][-1]), float(cs_exp_cb[count_cb][-1]))
                    ocs_list.append(line)
                    count_ca += 1
                    count_cb += 1
                if resn_ca > resn_cb:
                    line = '%4s %3s  %6.2f  %6.2f\n' % (cs_exp_cb[count_cb][0], cs_exp_cb[count_cb][1], 999.00, float(cs_exp_cb[count_cb][-1]))
                    ocs_list.append(line)
                    count_cb += 1
                if resn_ca < resn_cb:
                    line = '%4s %3s  %6.2f  %6.2f\n' % (cs_exp_ca[count_ca][0], cs_exp_ca[count_ca][1], float(cs_exp_ca[count_ca][-1]), 999.00)
                    ocs_list.append(line)
                    count_ca += 1
            except:
                break
        res_old = int(ocs_list[0].split()[0])
        count0 = 0
        count1 = 0
        safe = 0
        fd = open('%s.ocs' % cs_exp_name, 'w')
        while count0 < len(ocs_list):
            safe += 1
            if safe > len(cs_exp_ca)*5:
                break
            res_new = int(ocs_list[count0].split()[0])
            if res_old + count1 == res_new:
                fd.write(ocs_list[count0])
                count0 += 1
                count1 += 1
            else:
                fd.write('%4s UNK  999.00  999.00\n' % (res_old + count1))
                count1 += 1
        fd.close()
    except:
        fd = open('%s.ocs' % cs_exp_name, 'w')
        cs_file = open('%s' % cs_exp).readlines()
        reference = cs_file[0]
        for line in cs_file[1:]:
            fd.write(line)
        fd.close()
    return reference


def check_seq(residues, cs_exp_name):
    """ Compares if the aminoacidic sequence in the bmrb file matchs 
    the one in the pdb file."""
#read a pdb and extract the sequence using three leter code and save it to a list
    ok = 1
    ocslist = [] #contains the sequence from the ocs file
    ocslist_full = [] #contains the whole ocs file
    ocslist_full_new = [] #contains the corrected ocs file i.e. including UNK
    for line in open('%s.ocs' % (cs_exp_name)).readlines():
        ocslist_full.append(line)
        ocslist.append(line.split()[1])
    indelfirst, indellast = align(ocslist, residues)
    if indelfirst == 0 and indellast == 0:
        ocslist_full_new = list(ocslist_full)
    else:
        firstocs = int(ocslist_full[0].split()[0])
        lastocs = int(ocslist_full[-1].split()[0])
        newfirst = firstocs - indelfirst
        start = 0
        stop = len(ocslist_full)
        if indelfirst < 0:
            start = abs(indelfirst)
        if indellast < 0:
            stop  = len(ocslist_full) + indellast
        line = ('%s' % ocslist_full[0])
        for i in range(newfirst, firstocs): #works only if indelfirst is greater than 0
            line = ('%4s %3s  %6.2f  %6.2f\n' % (i, 'UNK', 999.00, 999.00))
            ocslist_full_new.append(line)
        for i in range(start, stop):
            line = ('%s' % ocslist_full[i])
            ocslist_full_new.append(line)
        for i in range(lastocs, lastocs+indellast):#works only if indellast is positive
            line = ('%4s %3s  %6.2f  %6.2f\n' % (i, 'UNK', 999.00, 999.00))
            ocslist_full_new.append(line)
#check if both sequences match
    fd = open('%s.ocs' % (cs_exp_name), 'w')
    for i, residue in enumerate(residues):
        if residue == ocslist_full_new[i].split()[1] or ocslist_full_new[i].split()[1] == 'UNK':
            fd.write('%s' % ocslist_full_new[i])
        else:
            pdb_res = residue
            ocs_res = ocslist[i-indelfirst]
            pdb_num = i + 1
            ocs_num = i + 1 - indelfirst
            ocs_seq = ' '.join(ocslist)
            ok = 0
            for ocs in ocslist:
                if ocs_seq.endswith('UNK'):
                    ocs_seq = ocs_seq[:-4]
                else:
                    break
            print 'The residue %s-%s in your PDB file does not match with residue %s-%s in the chemical shift file %s %s %s' % (pdb_res, pdb_num, ocs_res, ocs_num, ocslist[i-1], ocslist[i], ocslist[i+1])
            break
    fd.close()
    return ok, ocslist_full_new


def align(three0_list, three1_list):
    """This function aligns two sequences using a brute force algorithm. 
Returns how many positions the first sequence is shifted at the beginning 
and how many at the end. The sequences must have not indels and the first 
sequence must be shorter than the second"""

    three2one={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':
'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V', 'UNK':'U'}
#convert both lists from 3 letter code to one letter code
    one0_list = []
    one1_list = []
    for res in three0_list:
        one0_list.append(three2one[res])
    for res in three1_list:
        one1_list.append(three2one[res])
#convert one letter code list to strings
    a_seq = ''.join(one0_list)
    b_seq = ''.join(one1_list)
# get the length of the sequences
    len_a_seq = len(a_seq)
    len_b_seq = len(b_seq)
#create two new sequences of the same length
    seq_1 = len(b_seq) * '.' + a_seq
    seq_2 = b_seq + len(a_seq) * '.'
# compare both sequences, the trick is to delete (iteratively) the first chacter
# of sequence 1 and the last of sequence 2.
    dif = []
    for shift in range(0, len(seq_1)-1):
        seq_1 = seq_1[1:]
        seq_2 = seq_2[0:len(seq_2)-1]
        matching = 0
        for i in range(0, len(seq_1)):
            if seq_1[i] == seq_2[i]:
                matching += 1
        dif.append(matching)
    maximun = max(dif)
    for values in range(0, len(dif)):
        if maximun == dif[values]:
            beginning = len_b_seq-(values+1)
            endding = len_b_seq-len_a_seq-beginning
    return beginning, endding



def pose_from_pdb(pdb_file_name):
    """Gets information from the pdb like the number of residues, the sequence,
    the number of states and the name of the object"""
    pose = pdb_file_name
    remStr = "all and not (alt ''+A)"
    cmd.remove(remStr)
    cmd.alter('all', "alt=''")
    stored.residues = []
    stored.ResiduesNumber = []
    cmd.iterate('(name ca)','stored.residues.append((resn))')
    cmd.iterate('all','stored.ResiduesNumber.append((resi))')
    first = int(stored.ResiduesNumber[0])
    cmd.alter(pose, 'resi=str(int(resi)-%s)' % (first))
    cmd.sort(pose)
    states = cmd.count_states('all') + 1
    return pose, stored.residues, len(stored.residues), states



def get_phi(res_num, state):
    """Computes the phi torsional angle"""
    if res_num != 0:
        cmd.select('A', 'resi %s and name C' % (res_num-1))
        cmd.select('B', 'resi %s and name N' % res_num)
        cmd.select('C', 'resi %s and name CA' % res_num)
        cmd.select('D', 'resi %s and name C' % res_num)
        return cmd.get_dihedral('A', 'B', 'C', 'D', state)
    else:
        return float('nan')


def get_psi(res_num, state, total_residues):
    """Computes the psi torsional angle"""
    if res_num != total_residues - 1:
        cmd.select('A', 'resi %s and name N' % res_num)
        cmd.select('B', 'resi %s and name CA' % res_num)
        cmd.select('C', 'resi %s and name C' % res_num)
        cmd.select('D', 'resi %s and name N' % (res_num+1))
        psi = cmd.get_dihedral('A', 'B', 'C', 'D', state)
        return psi
    else:
        return float('nan')


def get_omega(res_num, state, total_residues):
    """Computes the omega torsional angle"""
    if res_num != total_residues-1:
        cmd.select('A', 'resi %s and name CA' % res_num)
        cmd.select('B', 'resi %s and name C' % res_num)
        cmd.select('C', 'resi %s and name N' % (res_num+1))
        cmd.select('D', 'resi %s and name CA' % (res_num+1))
        omega = cmd.get_dihedral('A', 'B', 'C', 'D', state)
        return omega
    else:
        return float('nan')


def get_chi1(res_num, res_name, state):
    """Computes the chi1 torsional angle"""
    if res_name not in ['ALA', 'GLY', 'PRO']:
        cmd.select('A', 'resi %s and name N' % res_num)
        cmd.select('B', 'resi %s and name CA' % res_num)
        cmd.select('C', 'resi %s and name CB' % res_num)
        cmd.select('D', 'resi %s and (name CG or name CG1 or name OG1 or name OG or name SG)' % res_num)
        chi1 = cmd.get_dihedral('A', 'B', 'C', 'D', state)
        return chi1
    else:
        return float('nan')


def get_chi2(res_num, res_name, state):
    """Computes the chi2 torsional angle"""
    if res_name not in ['ALA', 'GLY', 'PRO', 'SER', 'THR', 'VAL', 'CYS']:
        cmd.select('A', 'resi %s and name CA' % res_num)
        cmd.select('B', 'resi %s and name CB' % res_num)
        cmd.select('C', 'resi %s and (name CG or name CG1 or name OG1 or name OG)' % res_num)
        cmd.select('D', 'resi %s and (name CD or name CD1 or name OD1 or name ND1 or name SD)' % res_num)
        chi2 = cmd.get_dihedral('A', 'B', 'C', 'D', state)
        return chi2
    else:
        return float('nan')


def load(path):
    """Load the files containing the theoretical chemical shifts. Creates a 
    dictionary to store the data."""
    aminoacids = ['ALA','ARG','ASN','ASP','GLU','GLN','GLY','HIS','ILE','LEU',
    'LYS','MET','PHE','PRO','SER','THR','TYR','TRP','VAL']
    Db = {}
    for aminoacid in aminoacids:
        vector = []
        for line in open(os.path.join(path, 'CS_DB', 'CS_db_%s' % aminoacid)).readlines():
            vector.append(map(float, line.split()))
        Db[aminoacid] = vector
    return Db


def near_pro(omega, psi, Db):
    """Computes the chemical shifts from the psi and omega torsional angles
    by linear interpolation from the theoretical values stored in Db. 
    This funcion works only for proline"""
    points = []
    values_Ca = []
    values_Cb = []
    if omega  <= -90: # torsional angles are circular, for example -180=180. angles smaller than -90
        omega = 180   # are closer to 180 than to 0
    lista = np.array([0., 180.])# cheshift databse has only two values for proline omega angle, this
    index = (np.abs(lista-omega)).argmin() # two lines calculate the nearest omega angle, in the datbase,
    nearestOmega = lista[index] 
    for line in Db['PRO']: # PRO database is small. hence to calculate the theoretical CS i just take
        if line[0] == nearestOmega: # all the values with the nearestOmega
            vi, yi, zi = line[1], line[4], line[5]
            points.append(vi)
            values_Ca.append(yi), values_Cb.append(zi)
    points = np.array(points)
    values_Ca = np.array(values_Ca)
    values_Cb = np.array(values_Cb)
    values_Ca_New = griddata(points, values_Ca, (psi), method='linear') #linear interpolation
    values_Cb_New = griddata(points, values_Cb, (psi), method='linear') #linear interpolation
    return values_Ca_New, values_Cb_New


def near(phi, psi, chi1, chi2, res_name, Db):
    """Computes the chemical shifts from the torsional angles by linear 
    interpolation from the theoretical values stored in Db. 
    This funcion works for non-proline residues"""
    points = []
    values_Ca = []
    values_Cb = []
    phi_list = []
    phi_round = int(round(phi, -1)) # round to the nearest values in the database
    psi_round = int(round(psi, -1))
    chi1_rotamers =  np.array([-180., -150., -120., -90., -60., -30., 0., 30., 60., 90., 120., 150., 180.])
    index = (np.abs(chi1_rotamers-chi1)).argmin()
    nearestChi1_A = chi1_rotamers[index]
    chi1_rotamers_new = np.delete(chi1_rotamers, index)
    index = (np.abs(chi1_rotamers_new-chi1)).argmin()
    nearestChi1_B = chi1_rotamers_new[index]
    if phi > phi_round: # for phi and psi angles get the two nearest values in the database
        phi_range = range(phi_round, phi_round+20, 10)
    else:
        if phi_round == -180:
            phi_round = 180
        phi_range = range(phi_round-10, phi_round+10, 10)
    if psi > psi_round:
        psi_range = range(psi_round, psi_round+20, 10)
    else:
        if psi_round == -180:
            psi_round = 180
        psi_range = range(psi_round-10, psi_round+10, 10)
    for phi_value in phi_range: # A trick to avoid reading the whole list. The indexes where
        y = int(phi_value * 0.1 + 19)#  the necesarry values are stored (start and end) are calculated.
        if y > 37: 
            y = -(37-y)
        lenght = (len(Db[res_name])/37)
        end = lenght * y
        start = end-lenght
        for i in range(start, end):
            phi_list.append(Db[res_name][i])
    if res_name in ['ALA', 'GLY']:
        for line in phi_list:
            for psi_value in psi_range:
                if psi_value == line[1]:
                    ui, vi, yi, zi = line[0], line[1], line[4], line[5] 
                    vector = ui, vi
                    points.append(vector)
                    values_Ca.append(yi), values_Cb.append(zi)
        points = np.array(points)
        values_Ca = np.array(values_Ca)
        values_Cb = np.array(values_Cb)
        values_Ca_New = griddata(points, values_Ca, (phi, psi), method='linear')
        values_Cb_New = griddata(points, values_Cb, (phi, psi), method='linear')
        return  values_Ca_New, values_Cb_New
    elif res_name in ['SER', 'THR', 'VAL', 'CYS']:
        for line in phi_list:
            for psi_value in psi_range:
                if psi_value == line[1] and (line[2] == nearestChi1_A or line[2] == nearestChi1_B):
                    ui, vi, wi, yi, zi = line[0], line[1], line[2], line[4], line[5] 
                    vector = ui, vi, wi
                    points.append(vector)
                    values_Ca.append(yi), values_Cb.append(zi)
        points = np.array(points)
        points = np.array(points)
        values_Ca = np.array(values_Ca)
        values_Cb = np.array(values_Cb)
        values_Ca_New = griddata(points, values_Ca, (phi, psi, chi1), method='linear')
        values_Cb_New = griddata(points, values_Cb, (phi, psi, chi1), method='linear')
        return values_Ca_New, values_Cb_New
    else:
        lista = []
        for i in range (0, 3):
            rotamer = phi_list[i][3]
            if rotamer < 0:
                rotamer = rotamer + 360
            lista.append(rotamer)
        if 0. in lista:
            lista.append(360)
        lista = np.array(lista)
        if chi2 < 0:
            chi2 = chi2 + 360
        index = (np.abs(lista-chi2)).argmin()
        nearestChi2 = lista[index] 
        if nearestChi2 > 180:
            nearestChi2 = nearestChi2 - 360
        for line in phi_list:
            for psi_value in psi_range:
                if psi_value == line[1] and line[3] == nearestChi2 and (line[2] == nearestChi1_A or line[2] == nearestChi1_B):
                    ui, vi, wi, yi, zi = line[0], line[1], line[2], line[4], line[5] 
                    vector = ui, vi, wi
                    points.append(vector)
                    values_Ca.append(yi), values_Cb.append(zi)
        points = np.array(points)
        values_Ca = np.array(values_Ca)
        values_Cb = np.array(values_Cb)
        values_Ca_New = griddata(points, values_Ca, (phi, psi, chi1), method='linear')
        values_Cb_New = griddata(points, values_Cb, (phi, psi, chi1), method='linear')
        return values_Ca_New, values_Cb_New


def get_inout(state, residues, total_residues):
    """Computes if the phi-psi torsional angles of a residue belong to high or 
    low prababilities areas in the ramachandra plot. 
    Uses information derived from the Neighbor-dependent Ramachandran 
    Distributions http://dunbrack.fccc.edu/ndrd/ from Dunbrack Lab"""
    def myround(x, base=5):
        return int(base * round(x/base))

    inout_list = []
    for res_num in range(0, total_residues):
        if res_num == 0 or res_num == len(residues)-1:
            a = [residues[res_num], res_num, 'nan']
            inout_list.append(a)
        else:
            triple = residues[res_num-1]+residues[res_num]+residues[res_num+1]
            phi = get_phi(res_num, state)
            psi = get_psi(res_num, state, total_residues)
            phi_round = myround(phi)
            psi_round = myround(psi)
            fd = open(os.path.join(path, 'TRIPLE_5', '%s.dat' % triple)).readlines()
            if '%4s %4s\n' % (phi_round, psi_round) not in fd:
                a = [residues[res_num], res_num, 'low']
                inout_list.append(a)
            else:
                a = [residues[res_num], res_num, 'high']
                inout_list.append(a)
    return inout_list


def get_outliers(cs_exp_name): 
    """Computes if an observed chemical shift is within 3 standard deviation 
    from the mean of observed values. The observed values were taken from the 
    BMRB http://www.bmrb.wisc.edu/ref_info/statsel.htm"""
    boundaries_CA = {'ALA':(47.23,59.11), 'ARG':(49.86,63.72), 'ASP':(48.57,60.81), 'ASN':(47.88,59.22), 'CYS':(48.28,68.38), 'GLU':(51.07,63.61), 'GLN':(50.18,63.02), 'GLY':(41.37,49.35), 'HIS':(49.57,63.49), 'ILE':(53.58,69.72), 'LEU':(49.29,62.01), 'LYS':(50.38,63.58), 'MET':(49.42,62.86), 'PHE':(50.42,65.84), 'PRO':(58.79,67.91), 'SER':(52.50,64.98), 'THR':(54.47,70.01), 'TRP':(49.95,65.43), 'TYR':(50.57,65.75), 'VAL':(53.92,71.14)}
    boundaries_CB = {'ALA':(13.60,24.34), 'ARG':(25.17,36.15), 'ASP':(36.01,45.73), 'ASN':(33.71,43.67), 'CYS':(14.34,50.88), 'GLU':(24.88,35.08), 'GLN':(23.73,34.59), 'HIS':(24.01,36.43), 'ILE':(32.59,44.59), 'LEU':(36.70,47.86), 'LYS':(27.46,38.08), 'MET':(26.32,39.58), 'PHE':(33.74,46.16), 'PRO':(28.31,35.39), 'SER':(59.26,68.32), 'THR':(64.21,75.19), 'TRP':(23.89,36.07), 'TYR':(32.82,45.72), 'VAL':(27.37,38.05)}
    fd = open('%s.ocs' % (cs_exp_name))
    fd.readline()
    low_Ca = []
    high_Ca = []
    low_Cb = []
    high_Cb = []
    for line in fd:
        res_num = line.split()[0]
        res_name = line.split()[1]
        cs_exp_Ca_value = float(line.split()[2])
        cs_exp_Cb_value = float(line.split()[3])
        if cs_exp_Ca_value < 999.00:
            res_outlier = res_name.title() + res_num
            if  cs_exp_Ca_value < boundaries_CA[res_name][0]:
                low_Ca.append(res_outlier)
            elif  cs_exp_Ca_value > boundaries_CA[res_name][1]:
                high_Ca.append(res_outlier)
        if cs_exp_Cb_value < 999.00:
            res_outlier = res_name.title() + res_num
            if  cs_exp_Cb_value < boundaries_CB[res_name][0] :
                low_Cb.append(res_outlier)
            elif  cs_exp_Cb_value > boundaries_CB[res_name][1]:
                high_Cb.append(res_outlier)
    return low_Ca, high_Ca, low_Cb, high_Cb


def get_chemical_shifts(ocs_file, residues, total_residues, state, reference, Db):
    """Call the near and near_pro function only for the residues with observed 
    chemical shifts. Returns the list of computed chemical shifts and the list 
    of phi, psi, chi1 and chi2 torsional angles"""
    chemical_shifts = []
    res_num = 0
    tors_list = []
    for line in open(ocs_file).readlines(): 
        if line.split()[1] == 'UNK':
            a = ['UNK', 999.00, 999.00]
            chemical_shifts.append(a)
            res_name = residues[res_num]
            phi = get_phi(res_num, state)
            psi = get_psi(res_num, state, total_residues)
            chi1 = get_chi1(res_num, res_name, state)
            chi2  = get_chi2(res_num, res_name, state)
            res_name = residues[res_num]
            tors_list.append([res_name, phi, psi, chi1, chi2])
            res_num += 1
        else:
            res_name = residues[res_num]
            phi = get_phi(res_num, state)
            psi = get_psi(res_num, state, total_residues)
            chi1 = get_chi1(res_num, res_name, state)
            chi2  = get_chi2(res_num, res_name, state)
            try:
                res_name_next = residues[res_num+1]
            except:
                res_name_next = 'GLY'
            if res_name != 'PRO' and res_name != 'CYS':
                try:
                    values_Ca_New, values_Cb_New = near(phi, psi, chi1, chi2, res_name, Db)
                except:
                    values_Ca_New = 999.00
                    values_Cb_New = 999.00
            elif res_name == 'CYS':
                values_Ca_New = 999.00
                values_Cb_New = 999.00
            else:
                try:
                    omega = get_omega(res_num-1, state, total_residues)
                    values_Ca_New, values_Cb_New = near_pro(omega, psi, Db)
                except:
                    values_Ca_New, values_Cb_New = 999.00, 999.00 
            if res_name_next == 'PRO':
                a = [res_name, round((values_Ca_New -1.95 + reference), 2), round((values_Cb_New + reference),2)]
                chemical_shifts.append(a)
            else:
                a = [res_name, round((values_Ca_New + reference), 2), round((values_Cb_New + reference),2)]
                chemical_shifts.append(a)
            res_num += 1
            tors_list.append([res_name, phi, psi, chi1, chi2])
    return chemical_shifts, tors_list


def cs_2_colors(cs_exp_name, pose, residues, total_residues, states, reference, Db):
    """Calculates the theoretical chemical shifts and calculates the errors 
    |CS_theo-Cs_exp|. The errors are discretized and appended to the b-factor 
    column of the protein model. This function is exclusive of the validation 
    routines """
    cs_theo_list = []
    tors_matrix = []
    for state in range(1, states):
        cs_list, tors_list = get_chemical_shifts('%s.ocs' % (cs_exp_name), residues, total_residues, state, reference, Db)
        cs_theo_list.append(cs_list)
        tors_matrix.append(list(tors_list))
########## create list with experimental chemical shifts ###########
    exp_list = [[],[]]
    for nucleus in [0, 1]:
        for line in open('%s.ocs' % (cs_exp_name)).readlines():
            exp_list[nucleus].append(float(line.split()[nucleus+2]))
#############################################################
# for each conformation let find how to fix the side-chains #
#############################################################
    new_color_matrix = []
    for state in range(0, states-1):
        color_list = [[],[]]
        disc_list = [[],[]] 
        for nucleus in [0, 1]:
            if nucleus == 0:
                nucleus_name = 'Ca'
                cut = 1.45
            else:
                nucleus_name = 'Cb'
                cut = 1.77
            betalist = []
            betalist_disc = []
            for residue in range(0, len(exp_list[nucleus])):
                theo_value = cs_theo_list[state][residue][nucleus+1]
                if theo_value > 100:
                    betavalue = 999.00
                else:
                    betavalue = abs(theo_value - exp_list[nucleus][residue])
                betalist.append(betavalue)
            for value in betalist:
                if value > 100.:
                    value_disc = -2.0
                    color_list[nucleus].append('white')
                elif value >= (cut*2):
                    value_disc = -1.00
                    color_list[nucleus].append('red') 
                elif value >= cut:
                    value_disc = -0.50
                    color_list[nucleus].append('yellow')
                else:
                    value_disc = 0.00
                    color_list[nucleus].append('green')
                betalist_disc.append(value_disc)
                disc_list[nucleus] = betalist_disc
        inout_list = get_inout(state, residues, total_residues)
        kai2 = {'ARG':[-60,60,180],'ASN':[-75,-20,30],'ASP':[-15,0],'CYS':[float('nan')],'GLU':[-60,60,180],'GLN':[-65,65,180],'HIS':[-75,60,75],'ILE':[-60,60,180],'LEU':[65,175],'LYS':[-60,60,180],'MET':[-60,60,180],'PHE':[-90,90],'SER':[float('nan')],'THR':[float('nan')],'TYR':[-85,80],'TRP':[-105,-90,90],'VAL':[float('nan')]}
        report_green = []
        report_cyan = []
        report_all = []
        new_color = []
        for i in range(0, len(inout_list)):
            try:
                phi = tors_matrix[state][i][1]
                psi = tors_matrix[state][i][2]
                chi1 = tors_matrix[state][i][3]
                chi2 = tors_matrix[state][i][4]
            except:
                pass
            res_name = inout_list[i][0]
            res_num = inout_list[i][1]
            Ca_disc_value = disc_list[0][i]
            Cb_disc_value = disc_list[1][i]
            report_all.append('%5s%6s%9s%9s%11.4f%10.4f%10.4f%10.4f\n' % (res_name, res_num, color_list[0][i], color_list[1][i], phi, psi, chi1, chi2))
            if Ca_disc_value == 0.00 and Cb_disc_value == 0.00:
                new_color.append(0)
                if inout_list[i][2] == 'low':
                    report_green.append('%5s%6s%9s%9s%11.4f%10.4f%10.4f%10.4f\n' % (res_name, res_num, color_list[0][i], color_list[1][i], phi, psi, chi1, chi2))
            elif Ca_disc_value == -2.00 or Cb_disc_value == -2.00: # if something is missing check if everything is missing
                new_color.append(0)
            elif res_name in ['ALA', 'PRO']:
                new_color.append(0)
            else:
                if inout_list[i][2] == 'low': # does not make sense to try to fix this
                    new_color.append(0)
                else:
                    fix = 0
                    res_name_next = inout_list[i+1][0]
                    for chi1 in range(-180, 180, 30):
                        for chi2 in kai2[res_name]:
                            try:
                                values_Ca_New, values_Cb_New = near(phi, psi, chi1, chi2, res_name, Db)
                            except:
                                values_Ca_New = float('nan')
                                values_Cb_New = float('nan')
                            if res_name_next == 'PRO':
                                Theo_Ca = values_Ca_New - 1.95 + reference
                                Theo_Cb = values_Cb_New + reference
                            else:
                                Theo_Ca = values_Ca_New + reference
                                Theo_Cb = values_Cb_New + reference
                            if abs(Theo_Ca - exp_list[0][i]) <= 1.45 and abs(Theo_Cb - exp_list[1][i]) <= 1.77:
                                report_cyan.append('%5s%6s%9s%9s%11.4f%10.4f%10.4f%10.4f\n' % (res_name, res_num, color_list[0][i], color_list[1][i], phi, psi, chi1, chi2))
                                fix += 1
                    if fix == 0:
                        new_color.append(0)
                    else:
                        new_color.append(1)
            new_color_matrix.append(new_color)
    ############# write the ASCII report file #####################################
            fd = open('%s_%02d.sc' % (pose, state ), 'w')
            fd.write('*'*79+'\n')
            fd.write('HEADER    CheShift validation report of the protein %s, model %02d\n' % (pose, state))
            fd.write('*'*79+'\n')
            fd.write('REMARK The residues listed below, are those that occupy highly-populated\nREMARK regions of the Ramachandran map and for which a change in the side-chain\nREMARK torsional angles will lead to a good agreement between the observed and\nREMARK predicted chemical shifts for both the 13Ca and 13Cb nuclei, i.e., they\nREMARK will became "green" rather than yellow or red.\n')
            if len(report_cyan) != 0:
                fd.write('ResName Res# CA_color CB_color    phi       psi       chi1      chi2  \n')
                for i in report_cyan:
                        fd.write('%s' % i)
            else:
                fd.write('\n\n                      +-------------------------------+\n')
                fd.write('                      + There are no residues to show +\n')
                fd.write('                      +-------------------------------+\n\n\n')
            fd.write('*'*79+'\n')
            fd.write('REMARK The residues listed below belong to low-populated regions of the\nREMARK Ramachandran map and, hence, the good agreement in term of chemical\nREMARK shift (green) should be considered with caution because there is no\nREMARK one-to-one correspondence between value of the chemical shift and the\nREMARK conformation of the residue.\n')
            if len(report_green) != 0:
                fd.write('ResName Res# CA_color CB_color    phi       psi       chi1      chi2  \n')
                for i in report_green:
                        fd.write('%s' % i)
            else:
                fd.write('\n\n                      +-------------------------------+\n')
                fd.write('                      + There are no residues to show +\n')
                fd.write('                      +-------------------------------+\n\n\n')
            low_Ca, high_Ca, low_Cb, high_Cb = get_outliers(cs_exp_name) # get residues with unsual experimental chemical shifts.
            fd.write('*'*79+'\n')
            fd.write('REMARK The residues listed below have unusual observed chemical shifts (CS)\nREMARK values according to the statistics of the BMRB database.\nREMARK (http://www.bmrb.wisc.edu/ref_info/statsel.htm)\n')
            if len(low_Ca) > 0 or len(high_Ca) > 0 or len(low_Cb) > 0 or len(high_Cb) > 0:
                if len(low_Ca) > 0:
                    low_Ca = ", ".join(low_Ca)
                    fd.write('\nResidues with 13Ca CS below 3 standard deviation from the expected value\n%s\n' % low_Ca)
                if len(high_Ca) > 0:
                    high_Ca = ", ".join(high_Ca)
                    fd.write('\nResidues with 13Ca CS above 3 standard deviation from the expected value\n%s\n' % high_Ca)
                if len(low_Cb) > 0:
                    low_Cb = ", ".join(low_Cb)
                    fd.write('\nResidues with 13Cb CS below 3 standard deviation from the expected value\n%s\n' % low_Cb)
                if len(high_Cb) > 0:
                    high_Cb = ", ".join(high_Cb)
                    fd.write('\nResidues with 13Cb CS above 3 standard deviation from the expected value\n%s\n' % high_Cb)
            else:
                fd.write('\n\n                      +-------------------------------+\n')
                fd.write('                      + There are no residues to show +\n')
                fd.write('                      +-------------------------------+\n\n\n')
        
            fd.write('*'*79+'\n')
            fd.write('REMARK All residues are listed below\n')
            fd.write('ResName Res# CA_color CB_color    phi       psi       chi1      chi2  \n')
            for i in report_all:
                fd.write('%s' % i)
            fd.close()
    disc_list = [[],[]]
    color_list = [[],[]]
    for nucleus in [0,1]:
        if nucleus == 0:
            nucleus_name = 'Ca'
            cut = 1.45
        else:
            nucleus_name = 'Cb'
            cut = 1.77
        cs_theo_ave_list = []
        for residue in range(0, len(exp_list[nucleus])):
            i = 0
            lenght = 0
            for conformation in range(0, len(cs_theo_list)):
                if cs_theo_list[conformation][residue][nucleus+1] < 100.:
                    i = i + cs_theo_list[conformation][residue][nucleus+1]
                    lenght += 1
                else:
                    pass
            if lenght == 0:
                cs_theo_ave_list.append(-999.00)
            else:
                cs_theo_ave_list.append((i/lenght))
        betalist = []
        betalist_disc = []
        for i in range(0, len(exp_list[nucleus])):
            betavalue = abs(cs_theo_ave_list[i] - exp_list[nucleus][i])
            betalist.append(betavalue)
        for value in betalist:
            if value > 100:
                value_disc = -2.0
                color_list[nucleus].append('white')
            elif value >= (cut*2):
                value_disc = -1.00
                color_list[nucleus].append('red') 
            elif value >= cut:
                value_disc = -0.50
                color_list[nucleus].append('yellow')
            else:
                value_disc = 0.00
                color_list[nucleus].append('green')
            betalist_disc.append(value_disc)
            disc_list[nucleus] = betalist_disc
        for index in range(0, total_residues):
            cmd.alter('%s and resi %s' % (pose, index), 'b=%s' % betalist_disc[index])
        cmd.save('%s_%s.pdb' % (pose, nucleus_name), state=0)
    blue_list = []
    for column in range(0, len(new_color_matrix[0])):
        suma = 0
        for row in range(0 ,len(new_color_matrix)):
            suma += new_color_matrix[row][column]
        if suma > 0:
            blue_list.append(1)
        else:
            blue_list.append(0)
    new2_color = []
    for i in range(0, len(inout_list)):
        Ca_disc_value = disc_list[0][i]
        Cb_disc_value = disc_list[1][i]
        if blue_list[i] == 1:
            new2_color.append(2.00)
        elif Ca_disc_value == 0.00 and Cb_disc_value == 0.00:
            new2_color.append(0.00)
        elif Ca_disc_value == -2.00 or Cb_disc_value == -2.00: # if something is missing check if everything is missing
                if Ca_disc_value != -2.00:
                    new2_color.append(Ca_disc_value)
                elif Cb_disc_value != -2.00:
                    new2_color.append(Cb_disc_value)
                else:
                    new2_color.append(-2.00)
        else:
            if Ca_disc_value == -1.00 or Cb_disc_value == -1.00: 
                new2_color.append(-1.00)
            elif Ca_disc_value == -0.50 or Cb_disc_value == -0.50:
                new2_color.append(-0.50)
    for index in range(0, total_residues):
        cmd.alter('%s and resi %s' % (pose, index), 'b=%s' % new2_color[index])
    cmd.save('%s_CaCb.pdb' % pose, state=0)


def get_chemical_shifts_raw(residues, total_residues, state, Db):
    """Call the near and near_pro function for all the residues. This function 
    is exclusive of the prediction routines"""
    chemical_shifts = []
    for res_num in range(0, total_residues):
        res_name = residues[res_num]
        try:
            res_name_next = residues[res_num+1]
        except:
            res_name_next = 'GLY'
        if res_name != 'PRO' and res_name != 'CYS':
            try:
                phi = get_phi(res_num, state)
                psi = get_psi(res_num, state, total_residues)
                chi1 = get_chi1(res_num, res_name, state)
                chi2  = get_chi2(res_num, res_name, state)
                values_Ca_New, values_Cb_New = near(phi, psi, chi1, chi2, res_name, Db)
            except:
                values_Ca_New = 999.00
                values_Cb_New = 999.00
        elif res_name == 'CYS':
            values_Ca_New = 999.00
            values_Cb_New = 999.00
        else:
            try:
                omega = get_omega(res_num-1, state, total_residues)
                values_Ca_New, values_Cb_New = near_pro(omega, psi, Db)
            except:
                values_Ca_New, values_Cb_New = 999.00, 999.00 
        if res_name_next == 'PRO':
            a = [res_name, round((values_Ca_New -1.95), 2), round((values_Cb_New),2)]
            chemical_shifts.append(a)
        else:
           a = [res_name, round((values_Ca_New), 2), round((values_Cb_New),2)]
           chemical_shifts.append(a)
    return chemical_shifts


def raw(pose, residues, total_residues, states, Db):
    """Calculates the theoretical chemical shifts. This function is exclusive 
    of the prediction routine"""
    cs_theo_list = []
    for state in range(1, states):
        cs_list = get_chemical_shifts_raw(residues, total_residues, state, Db)
        cs_theo_list.append(cs_list)
    fd = open('%s.txt' % pose, 'w')
    fd.write('Ca Chemical Shifts\n')
    for residue in range(0, total_residues): 
        cs_theo_line = []
        for conformation in range(0, len(cs_theo_list)):
            res_name, Ca_shift, Cb_shift = cs_theo_list[conformation][residue]
            if float(Ca_shift) > 100.:
                Ca_shift = 999.00
            cs_theo_line.append('%6.2f' % (Ca_shift))
        res_line = "\t".join(cs_theo_line)
        fd.write('%s\t %s\n' % (residues[residue], res_line))
    fd.write('\nCb Chemical Shifts\n')
    for residue in range(0, len(residues)): 
        cs_theo_line=[]
        for conformation in range(0, len(cs_theo_list)):
            res_name, Ca_shift, Cb_shift = cs_theo_list[conformation][residue]
            if float(Cb_shift) > 100.:
                Cb_shift = 999.00
            cs_theo_line.append('%6.2f' % (Cb_shift))
        res_line = "\t".join(cs_theo_line)
        fd.write('%s\t %s\n' % (residues[residue], res_line))
    fd.close()


def clean(pose):
    """Deletes everything and load the validated models"""
    cmd.delete('all')
    cmd.load('%s_Ca.pdb' % pose)
    cmd.load('%s_Cb.pdb' % pose)
    cmd.load('%s_CaCb.pdb' % pose)
    cmd.intra_fit('%s_Ca' % pose, 0, quiet=1)
    cmd.intra_fit('%s_Cb' % pose, 0, quiet=1)
    cmd.intra_fit('%s_CaCb' % pose, 0, quiet=1)
    cmd.dss('all')
    cmd.disable('%s_Ca' % pose)
    cmd.disable('%s_Cb' % pose)

