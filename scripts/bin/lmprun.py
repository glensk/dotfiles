#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse

def _initialize_simulation(lmp):
    lmp.command("clear")
    lmp.command("units metal")
    lmp.command("dimension 3")
    lmp.command("boundary p p p")
    lmp.command("atom_style atomic")
    lmp.command("atom_modify map array")

def _create_lattice_command(lattice_const, lattice):
    lattice_cmd = "lattice custom {} ".format(lattice_const)
    lattice_cmd += "a1 {} {} {} a2 {} {} {} a3 {} {} {}".format(
        *lattice.flatten()
    )
    lattice_cmd += " basis 0.0 0.0 0.0"
    lattice_cmd += " basis 0.0 0.5 0.5"
    lattice_cmd += " basis 0.5 0.0 0.5"
    lattice_cmd += " basis 0.5 0.5 0.0"
    return lattice_cmd

def _create_region_command(lattice_const, lattice):
    lattice = lattice * lattice_const
    xhi = lattice[0,0]
    xy = lattice[1,0]
    yhi = lattice[1,1]
    xz = lattice[2,0]
    yz = lattice[2,1]
    zhi = lattice[2,2]

    region_cmd = "region box prism 0 {} 0 {} 0 {} ".format(
        xhi, yhi, zhi
    )
    region_cmd +=" {} {} {} units box".format(xy, xz,yz)
    return region_cmd

def _create_atoms(lmp, lattice_const, lattice):
    print('444')
    lattice_cmd = _create_lattice_command(lattice_const, lattice)
    print('555')
    print(lattice_cmd)
    lmp.command(lattice_cmd)
    print('666')

    region_cmd = _create_region_command(lattice_const, lattice)
    print('777')
    lmp.command(region_cmd)
    print('888')

    lmp.command("create_box 1 box")
    print('999')

    lmp.command(lattice_cmd)
    print('xxx')
    #lmp.command("mass 1 24.305")
    #lmp.command("mass 2 26.9815385")
    #lmp.command("mass 3 28.0855")
    #lmp.command("create_atoms 2 box")
    print('yyy')

def _define_interatomic_potential(lmp):
    #lmp.command("pair_style meam")
    #lmp.command("pair_coeff * * library.meam Al meam.alsimgcufe Al")

    lmp.command("pair_style eam/alloy")
    lmp.command("pair_coeff * * Al.eam.alloy_cutoff5_seed_914702_std_1.05861 Al")

    #print('uuu')
    #lmp.command("pair_style runner dir \"/Users/glensk/Dropbox/Albert/scripts/dotfiles/scripts/potentials/runner_v3ag_5000\" showewsum 1 showew yes resetew no maxew 1000000")
    #print('iii')
    #lmp.command("pair_coeff * * 7.937658735")
    #lmp.command("neighbor 1.0 nsq")
    #lmp.command("neigh_modify once no every 1 delay 0 check yes")

    lmp.command("neighbor 2.0 bin")
    lmp.command("delete_atoms overlap 0.3 all all")

def _define_settings(lmp):
    lmp.command("compute eng all pe/atom")
    lmp.command("compute eatoms all reduce sum c_eng")

def _run_minimization(lmp):
    lmp.command("reset_timestep 0")
    lmp.command("thermo_style custom step pe pxx pyy pzz pyz pxz pxy")
    lmp.command("dump 1 all atom 1 lammps.dump")
    lmp.command("min_style cg")
    lmp.command("minimize 1e-25 1e-25 5000 10000")

    #lmp.command("run 0")
def _print_compliance_components(lmp, C_1st_component):
    pressure_labels = ["pxx","pyy","pzz","pyz","pxz","pxy"]
    e = 0.002
    # debug
    #for x in pressure_labels:
    #    print "{} {}".format(x, lmp.get_thermo(x))

    for i in range(len(pressure_labels)):
        pressure_value = lmp.get_thermo(pressure_labels[i])
        barr2Gpa = 0.0001
        Compliance_value = -1*pressure_value/e * barr2Gpa
        C_2nd_component = i+1
        C_label = "{}{}".format(C_1st_component, C_2nd_component)
        print("{} : {} GPa".format(C_label, Compliance_value))

def run_fcc(lmp, lattice_const, lattice):
    print('111')
    _initialize_simulation(lmp)
    print('222')
    #lmp.command("read_data \"pos.lmp\"")
    _create_atoms(lmp, lattice_const, lattice)
    print('333')
    #sys.exit()

    _define_interatomic_potential(lmp)
    print('555')
    _define_settings(lmp)
    _run_minimization(lmp)
    return lmp

def _collect_strain_energies(lmp, strain_definition, strain_range):
    energies = []
    for strain_amount in strain_range:
        applied_strain = strain_definition * strain_amount
        lmp = _run_lammps_at_strain(lmp,
                                   e1=applied_strain[0],
                                   e2=applied_strain[1],
                                   e3=applied_strain[2],
                                   e4=applied_strain[3],
                                   e5=applied_strain[4],
                                   e6=applied_strain[5],
                                   )
        energy = lmp.get_thermo("pe")
        energies.append(energy)
    energies = np.array(energies)
    return energies

def _fitto_Ax2(x_vals, y_vals):
    A = np.linalg.lstsq(x_vals, y_vals, rcond=-1)[0]
    return A

def _find_compliance_viaenergy(lmp,
                               strain_definition,
                               strain_range,
                               V0):

    energies = _collect_strain_energies(lmp,
                                        strain_definition,
                                        strain_range)

    del_energies = np.vstack(energies - np.min(energies))
    strain_range_sqrd = np.vstack(strain_range ** 2)
    A = _fitto_Ax2(strain_range_sqrd, del_energies)

    evAng2GPa = 160.21766208
    compliance = A * 2*evAng2GPa/V0

    return compliance
