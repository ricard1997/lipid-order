"""
Unit and regression test for the lipid_order package.
"""

# Import package, test suite, and other packages as needed
#import lipid_order
#import pytest
import MDAnalysis as mda
import sys
from lipid_order import sn1, sn2
import matplotlib.pyplot as plt
def test_lipid_order_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "lipid_order" in sys.modules





def test_sn1():
    top = "./../../../../trajectories/0model/rep0/centered_prot.gro"
    traj = "./../../../../trajectories/0model/rep0/centered_prot.xtc"

    u = mda.Universe(top, traj)

    order_parameters = sn1(u, "(resname DSPC)", lipid = "DSPC", n_chain = 17, step = 100)

    plt.plot(order_parameters[0])
    plt.show()
    
def test_sn2():  
    top = "./../../../../trajectories/0model/rep0/centered_prot.gro"
    traj = "./../../../../trajectories/0model/rep0/centered_prot.xtc"  
    u = mda.Universe(top, traj)
    order_parameters = sn2(u, "(resname DSPC)", lipid = "DSPC", n_chain = 17, step = 100)

    plt.plot(order_parameters[0])
    plt.show()