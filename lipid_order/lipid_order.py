import MDAnalysis
import numpy as np

def get_vectors(lista):
    angles = []
    for i in range(len(lista)-1):
        vectores =lista[i+1].positions - lista[0].positions
        angles.append(vectores)
    angles = np.concatenate(angles, axis = 0)
    costheta = angles[:,2]**2/np.linalg.norm(angles, axis=1)**2
    order = np.mean(costheta)
    return order

def order_sn1(sel, lipid, n_chain):
    chains = []
    for i in range(n_chain):
        selections = [
                f"name C3{i+2}",
                f"name H{i+2}X and not name HX",
                f"name H{i+2}Y and not name HY",
                f"name H{i+2}Z and not name HZ",
        ]
        lista = []
        for selection in selections:
            atoms = sel.select_atoms(selection)
            if atoms.n_atoms != 0:
                lista.append(atoms)
        chains.append(get_vectors(lista))
    chains = np.array(chains)

    return chains


def order_sn2(sel, lipid, n_chain):
    chains = []
    for i in range(n_chain):
        selections = [
                f"name C2{i+2}",
                f"name H{i+2}R and not name HR",
                f"name H{i+2}S and not name HS",
                f"name H{i+2}T and not name HT",
        ]
        if lipid == "POPE" or lipid == "POPS":
            if selections[0] == "name C29":
                selections[1] == "name H91"
            if selections[0] == "name C210":
                selections[1] == "name H101"
        lista = []
        for selection in selections:
            atoms = sel.select_atoms(selection)
            if atoms.n_atoms != 0:
                lista.append(atoms)
        chains.append(get_vectors(lista))
    chains = np.array(chains)
    return chains

def sn1(u, sel_string, lipid, n_chain, start = 0, stop = -1, step = 1):
    order_parameters = []

    for ts in u.trajectory[start:stop:step]:
        sel = u.select_atoms(sel_string)
        order_parameters.append(order_sn1(sel, lipid, n_chain))
    order_parameters = np.array(order_parameters)
    serror = 1.5 * np.std(order_parameters, axis = 0)/np.sqrt(len(order_parameters))
    order_parameters = np.abs(1.5 * np.mean(order_parameters, axis = 0) - 1)
    return order_parameters, serror


def sn2(u, sel_string, lipid, n_chain, start = 0, stop = -1, step = 1):
    order_parameters = []

    for ts in u.trajectory[start:stop:step]:
        sel = u.select_atoms(sel_string)
        order_parameters.append(order_sn2(sel, lipid, n_chain))
    order_parameters = np.array(order_parameters)
    serror = 1.5 * np.std(order_parameters, axis = 0)/np.sqrt(len(order_parameters))
    order_parameters = np.abs(1.5 * np.mean(order_parameters, axis = 0) - 1)
    return order_parameters, serror


