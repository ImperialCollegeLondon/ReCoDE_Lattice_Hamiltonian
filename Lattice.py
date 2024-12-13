# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 15:53:16 2024

@author: ljh3218
"""
import scipy.constants as sp
import numpy as np
import pandas as pd
import math

import Redfield as Redfield
import Recombination as Recombination

from scipy.sparse import csr_array
from scipy import linalg
from copy import deepcopy
from itertools import combinations
from itertools import product
from dataclasses import dataclass, field

@dataclass
class Parameters:
    """Stores global parameters of the system"""
    
    temp: float
    """The temperature of the lattice in Kelvin."""
    e_peak: float
    """The energy of the peak of the spectral density function in eV. Typically 0.16 eV for organic molecules."""
    lambda_outer: float
    """The outer reorganisation energy of each molecule in the lattice in eV."""
    lambda_inner:float
    """The inner reorganisation energy of each molecule in the lattice in eV."""
    j0: float
    """A parameter used in determining the shape of the Mataga potenial. Units are eV."""
    r0j: float
    """A parameter used in determining the shape of the Mataga potenial. Units are angstrom."""
    e_singlet: float
    """The energy of a singlet exciton in eV. This should be lower than the bandgap of the material."""
    const_recombination: bool 
    """A switch which determines whether recombination rates are treated as a constant (True) 
    or are calculated using generalised Marcus-Levich-Jortner (False)."""

    @property 
    def kT(self):
        return sp.k/sp.e * self.temp

@dataclass
class Site:
    """Stores the properties of each lattice site"""
    
    coordinate: list[float]
    """A list of the (x,y,z) coordinates of the lattice site. Units are angstrom. """
    HOMO: float
    """Energy of the lattice site's highest occupied molecular orbital in eV."""
    LUMO: float
    """Energy of the lattice site's lowest unoccupied molecular orbital in eV."""
    id: int
    """The ID of the lattice site."""
    const_recombination: bool
    """A switch which determines whether recombination rates are treated as a constant (True) 
    or are calculated using generalised Marcus-Levich-Jortner (False)."""
    transition_dipole_ex: float
    """Number between zero and one which controls the rate at which excitons are generated in the lattice."""
    transition_dipole_ct: float 
    """Number between zero and one which controls the rate at which nearest neighbour 
    charge transfer states are generated in the lattice."""
    nearest_neighbour: list = field(default_factory = list)
    """List containing the IDs of all the lattice sites to include when calculating the size of the dipole-dipole coupling."""
    LUMO_coupling: list = field(default_factory = list)
    """List containing the electronic coupling between the LUMO of the lattice site and the LUMOs of all its adjacent neighbours."""
    HOMO_coupling: list = field(default_factory = list)
    """List containing the electronic coupling between the HOMO of the lattice site and the HOMOs of all its adjacent neighbours."""
    dipole_coupling: list = field(default_factory = list)
    """List containing the dipole coupling between the current lattice site and all 
    the lattice sites included in the list Nearest_Neighbour."""
    krec_ex: float = field(default=0, kw_only=True)
    """The decay rate of excitonic states in s-1."""
    krec_ct: float = field(default=0, kw_only=True)
    """The decay rate of nearest neighbour charge transfer states in s-1."""
    v_ct: float = field(default=0, kw_only=True)
    """The strength of the coupling between the ground state and the nearest neighbour charge transfer state in eV."""
    v_ex: float = field(default=0, kw_only=True)
    """The strength of the coupling between the ground state and the exciton state in eV."""

    # Check they are non-zero
    def __post_init__(self):
        if self.const_recombination: 
            if self.krec_ex == 0 or self.krec_ct == 0:
                raise ValueError("Recombination rates for exciton and CT states cannot be zero. Set krec_ex and krec_ct to finite values.")
        else:
            if self.v_ct == 0 or self.v_ex == 0:
                raise ValueError("Recombination rates for exciton and CT states cannot be zero. Set v_ex and v_ct to finite values.")

class Lattice:
    def __init__(self):
        self.sites=[]
        self.states=[]   
        
    def generate_uniform(self, size:int, HOMO:int, LUMO:int, 
                         dist_sites:int, min_dist_near_neighbour:float,
                         t0_homo:float, t0_lumo:float, 
                         d0:float, r0d:float, 
                         ct_abs_prop:float,
                         krec_ex:float = 0, krec_ct:float = 0, 
                         v_ex:float = 0, v_ct:float = 0,
                         const_recombination:bool = True):
        counter_id = 0
        self.const_recombination = const_recombination       
           
        site_vec = np.zeros((size**2, 3))
                
        for x,y in product(range(size), repeat = 2):
            #The 1 here is the transition dipole of the exciton
            if const_recombination:
                self.sites.append(Site(np.array([x*dist_sites, y*dist_sites, 0]),
                                       HOMO, LUMO, counter_id, const_recombination, 
                                       1, ct_abs_prop, 
                                       krec_ex = krec_ex, krec_ct = krec_ct))
            else:
                self.sites.append(Site(np.array([x*dist_sites, y*dist_sites, 0]),
                                       HOMO, LUMO, counter_id, const_recombination, 
                                       1, ct_abs_prop, 
                                       v_ct = v_ct, v_ex = v_ex))
            site_vec[counter_id,:] = np.array([x*dist_sites, y*dist_sites, 0])
            counter_id = counter_id + 1  
        self.site_vec = site_vec
        for ii,jj in combinations(range(counter_id),2):
            #Find distance from site ii to site jj
            distance_bet_sites = math.dist(self.sites[ii].coordinate, self.sites[jj].coordinate)
            #Only couple sites which are closer than whatever the nearest neighbour distance is 
            if distance_bet_sites < min_dist_near_neighbour:
                self.sites[ii].nearest_neighbour.append(self.sites[jj].id)
                self.sites[jj].nearest_neighbour.append(self.sites[ii].id)
                if math.isclose(distance_bet_sites-dist_sites, 0, abs_tol = 1e-3):
                    self.sites[ii].LUMO_coupling.append(t0_lumo)
                    self.sites[jj].LUMO_coupling.append(t0_lumo)
                    self.sites[jj].HOMO_coupling.append(t0_homo)
                    self.sites[ii].HOMO_coupling.append(t0_homo)
                else:
                    self.sites[ii].LUMO_coupling.append(0)
                    self.sites[jj].LUMO_coupling.append(0)
                    self.sites[jj].HOMO_coupling.append(0)
                    self.sites[ii].HOMO_coupling.append(0)
                self.sites[jj].dipole_coupling.append(d0*((distance_bet_sites-dist_sites)/r0d+1)**-3)
                self.sites[ii].dipole_coupling.append(d0*((distance_bet_sites-dist_sites)/r0d+1)**-3) 

    def sites_to_dataframe(self, nonuniform = 0):
        return pd.DataFrame.from_records([s.to_dict(nonuniform) for s in self.sites])
    
    def build_ham(self, params, F, min_dist_near_neighbour, dist_cs_min, disorder_site_ene, random_seed = 0):
        row = []
        col = []
        data = []
        basis = []
        const_recombination = self.const_recombination
        if random_seed != 0:
            np.random.seed(seed = random_seed)
        #Create empty lists to store the values of various quantities of interest 
        transdip_vec_ex = []
        krec_vec_ex = []
        transdip_vec_ct = []
        krec_vec_ct = []
        dist_he = []
        e_singlet = params.e_singlet
        is_ex = []    
        counter = 0
        for elec_pos in self.sites:
            for hol_pos in self.sites: 
                distance_ele_hole = math.dist(elec_pos.coordinate,hol_pos.coordinate)
                dist_vec = (elec_pos.coordinate - hol_pos.coordinate)
                dist_he.append(distance_ele_hole)
                coul_j = params.j0/(1+distance_ele_hole/params.r0j)
                #Setting generation and decay rates for CT and exciton states
                if distance_ele_hole == 0:
                    transdip_vec_ex.append(elec_pos.transition_dipole_ex)
                    is_ex.append(1)
                    if const_recombination:
                        krec_vec_ex.append(elec_pos.krec_ex)
                    else:
                        krec_vec_ex.append(elec_pos.v_ex)                     
                else:
                    transdip_vec_ex.append(0)
                    is_ex.append(0)
                    krec_vec_ex.append(0)
                if hol_pos.id in elec_pos.nearest_neighbour:
                    transdip_vec_ct.append(elec_pos.transition_dipole_ct)
                    #This is making the assumption that the coupling for CT/polaron states which aren't nearest neighbours
                    #to the ground is negligible. 
                    if const_recombination:
                        krec_vec_ct.append(elec_pos.krec_ct)
                    else: 
                        krec_vec_ct.append(elec_pos.v_ct)
                else:
                    transdip_vec_ct.append(0)
                    krec_vec_ct.append(0)                    
                #Calculate the values of the diagonal terms which go into the matrix  
                row.append(counter)
                col.append(counter)
                #Static disorder term added on as a random offset in the bandgap energy
                if distance_ele_hole == 0:
                    data.append(e_singlet + disorder_site_ene*np.random.standard_normal(1)[0])                     
                else:
                    data.append(elec_pos.LUMO - hol_pos.HOMO + disorder_site_ene*np.random.standard_normal(1)[0] - coul_j + np.dot(dist_vec, F))
                counter += 1
                basis.append([elec_pos.id, hol_pos.id])
        
        #off diagonal elements 
        #using combinations so will never have case where ii = jj
        for ii,jj in combinations(range(len(basis)),2):
            #CT state/exciton on basis[ii] and basis[jj] (exciton can only be on one of the two sites)
            if basis[jj][0] in self.sites[basis[ii][0]].nearest_neighbour and basis[ii][1] == basis[jj][1]:
                #index() returns the index of the first occurance of the specified value 
                counter_neighbour_ele = self.sites[basis[ii][0]].nearest_neighbour.index(basis[jj][0])                        
                row.append(ii)
                col.append(jj)
                #hopping for electron 
                data.append(self.sites[basis[ii][0]].LUMO_coupling[counter_neighbour_ele])
                #same again as H is a symmetric matrix
                row.append(jj)
                col.append(ii)
                data.append(self.sites[basis[ii][0]].LUMO_coupling[counter_neighbour_ele])
            #two excitons on sites which are close enough together that their interaction is non-zero
            #only exciton-exciton interaction considered currently is dipole-dipole
            if basis[ii][0] == basis[ii][1] and basis[jj][0] == basis[jj][1] and basis[jj][0] in self.sites[basis[ii][0]].nearest_neighbour:
                counter_neighbour_ele = self.sites[basis[ii][1]].nearest_neighbour.index(basis[jj][1])
                row.append(ii)
                col.append(jj)
                data.append(self.sites[basis[ii][0]].dipole_coupling[counter_neighbour_ele])
                row.append(jj)
                col.append(ii)
                data.append(self.sites[basis[ii][0]].dipole_coupling[counter_neighbour_ele])
            #CT state/exciton on basis[ii] and basis[jj] (exciton can only be on one of the two sites) 
            if basis[jj][1] in self.sites[basis[ii][1]].nearest_neighbour and basis[ii][0]==basis[jj][0]:  
                counter_neighbour_hole = self.sites[basis[ii][1]].nearest_neighbour.index(basis[jj][1])                      
                row.append(ii)
                col.append(jj)
                data.append(self.sites[basis[ii][1]].HOMO_coupling[counter_neighbour_hole])
                row.append(jj)
                col.append(ii)
                data.append(self.sites[basis[ii][1]].HOMO_coupling[counter_neighbour_hole])
        self.ham = csr_array((data, (row, col)), shape=(len(self.sites)**2,len(self.sites)**2))
        #add more properties to the lattice to be used in rate calculations 
        self.basis = basis
        self.transdip_vec_ct = transdip_vec_ct
        self.transdip_vec_ex = transdip_vec_ex
        self.krec_vec_ex = krec_vec_ex
        self.krec_vec_ct = krec_vec_ct         
        self.dist_he = dist_he
        self.is_ex = is_ex

    #Find the eigenvalues and eigenvectors of the Hamiltonian 
    def states_from_ham(self, params, max_energy_diff):  
        lambda_inner, lambda_outer = params.lambda_inner, params.lambda_outer, 
        evals,evecs = linalg.eigh(self.ham.toarray())
        #evals,evecs = eigsh(self.Ham, k=(len(self.sites)**2)-1, which = 'SM')
        evals = np.real(evals)
        evecs = np.real(evecs)
        #remove eigenvectors corresponding to states above the energy cut off
        evecs = evecs[:,evals<evals.min() + max_energy_diff]
        #remove eigenvalues corresponding to states above the energy cut off
        evals = evals[evals<evals.min() + max_energy_diff]
        #Declare a whole ton of empty lists which will be populated by the function 
        list_evecs = []
        dis_st = []
        IPR = []
        ex_char = [] 
        transdip_ex = [] 
        transdip_ct = [] 
        occupation_prob = [] 
        krec_ex = [] 
        krec_ct = []
        is_ex = np.array(self.is_ex)
        basis = np.array(self.basis)
        if not params.const_recombination:
            v_eff_ex = [] 
            v_eff_ct = []
        for i in range(len(evals)):
            list_evecs.append(evecs[:,i])
            occupation_probability = evecs[:,i]**2
            occupation_prob.append(occupation_probability)
            dis_st.append(occupation_probability @ self.dist_he)
            transdip_ex.append(occupation_probability @ self.transdip_vec_ex)
            transdip_ct.append(occupation_probability @ self.transdip_vec_ct)
            ex_char.append(occupation_probability @ is_ex)
            if params.const_recombination:
                IPR.append(calc_IPR(evecs[:,i], basis, is_ex))
                krec_ex.append(occupation_probability @ self.krec_vec_ex)
                krec_ct.append(occupation_probability @ self.krec_vec_ct)
            #Using generalised MLJ following Taylor and Kassal 2018/ D'Avino et al. J. Phys. Chem. Lett. 2016, 7, 536-540
            else:
                IPR.append(calc_IPR(evecs[:,i], basis, is_ex))
                effective_coupling_ex = evecs[:,i]@self.krec_vec_ex
                v_eff_ex.append(effective_coupling_ex)
                effective_inner_lambda_ex = lambda_inner*(1/IPR[i][0])
                effective_outer_lambda_ex = lambda_outer*(1/IPR[i][0])
                krec_ex.append(Recombination.decay_rate(effective_inner_lambda_ex, params.e_peak, effective_outer_lambda_ex, 
                                                         evals[i], effective_coupling_ex))

                effective_coupling_ct = evecs[:,i]@self.krec_vec_ct
                v_eff_ct.append(effective_coupling_ct)
                effective_inner_lambda_ct = lambda_inner*(1/IPR[i][1] + 1/IPR[i][2])
                effective_outer_lambda_ct = lambda_outer*(1/IPR[i][1] + 1/IPR[i][2])
                krec_ct.append(Recombination.decay_rate(effective_inner_lambda_ct, params.e_peak, effective_outer_lambda_ct, 
                                                         evals[i], effective_coupling_ct))                                
        states = pd.DataFrame({'energies': evals,'dis_eh': dis_st, 
                                'IPR': IPR, 'ex_char': ex_char,
                                'transdip_ex': transdip_ex,'transdip_ct': transdip_ct,
                                'krec_ex': krec_ex,'krec_ct': krec_ct,
                                'occupation_probability': occupation_prob})          
        #make it so the states are numbered 1 to n, not 0 to n-1
        states = states.sort_values(by=['energies'])
        states.insert(0, "state", np.arange(1, len(evals)+1))
        if params.const_recombination:
            states = states.assign(v_effective_ex = v_eff_ex, v_effective_ct = v_eff_ct)
        if len(self.states) == 0:
            self.states = states   
        else:
            #adding newly calculated values onto the pre-existing states dataframe 
            self.states = pd.concat([self.states, states], ignore_index=True)  
        return np.asarray(list_evecs)
    
    def get_rate_mat(self, params):
        '''Redfield rates using a secular approximation as in Marcus and Renger 2002 & Quantum biology revisited 2020'''
        e_peak, kT = params.e_peak, params.kT
        lambda_outer, lambda_inner, = params.lambda_outer, params.lambda_inner
        lambda_inner = Redfield.norm_J_Renger(lambda_outer+lambda_inner, e_peak)
        energies = np.array(self.states.energies)
        num_states = len(energies)
        occupation_prob = np.array([self.states.occupation_probability[i] for i in range(num_states)])        
        rates_mat = np.zeros((num_states, num_states), dtype = np.float32)
        inds = np.tril_indices_from(rates_mat, k=-1)
        len_inds = int(((num_states)*(num_states-1))/2) 
        w = np.array([energies[inds[0][k]] - energies[inds[1][k]] for k in range(len_inds)])
        C = Redfield.C_re_2D_array(w, lambda_inner, e_peak, kT) 
        #Here I pick the kth element from inds[0], which gives the row index, and the kth element from inds[1], which gives the column index
        #This corresponds to a transition between the inds[0][k] and inds[1][k] eigenstates
        gamma_Cw_real_tot = np.array([occupation_prob[inds[0][k]]@occupation_prob[inds[1][k]]*C[k] for k in range(len_inds)])
        rates_mat[inds] = gamma_Cw_real_tot
        rates_mat[(inds[1], inds[0])] = rates_mat[inds]*np.exp(-w/params.kT)  
        self.rates_mat = rates_mat

    def solve_steady(self,params):
        isteady_n = np.zeros(len(self.states))
        a = deepcopy(self.rates_mat)
        a = np.transpose(a)
        self.states['gen'] = self.states.transdip_ex + self.states.transdip_ct        
        for i in range(len(self.states)):
            a[i,i] = -sum(a[:,i])
            a[i,i] = a[i,i] - self.states.iloc[i].krec_ex - self.states.iloc[i].krec_ct
            isteady_n[i] = self.states.iloc[i].gen
        b = -isteady_n
        z = linalg.solve(a,b)
        self.states['z'] = z
    
#Code needs optimising 
def calc_IPR(eigvec, basis, is_ex):
        #These are calculated following method in D'Avino et al. J. Phys. Chem. Lett. 2016, 7, 536-540
        #Although I think there is a mistake in this paper - missing a squared sign???
        #See workbook 3, 14/11/23
        #Very slow currently as there are two for loops
        #Three if you count the for loop this is nested inside in get states
        num_states = len(eigvec)
        num_sites = int(basis[-1][0]+1)
        pop = eigvec**2
        
        frac_ex = sum(pop*is_ex)
        frac_ct = sum(pop*(1-is_ex))
        IPR_ex = frac_ex**2/sum((pop*is_ex)**2)

        #loop over lattice sites
        e_pop = np.zeros(num_sites)
        h_pop = np.zeros(num_sites)
        for k in range(num_sites):
            h_temp = 0
            e_temp = 0
            c_h_temp = 0
            c_e_temp = 0            
            #loop over basis elements
            for j in range(num_states):
                #if electron is on site k, add to sum as need to find probability electron on k for any hole position
                if basis[j][0] == k:
                    #Put this line here to avoid including basis elements which correspond to excitons in the sum 
                    if basis[j][0] != basis[j][1]:
                        e_temp += pop[j]
                    c_e_temp += pop[j]
                if basis[j][1] == k:
                    if basis[j][0] != basis[j][1]:
                        h_temp += pop[j]
                    c_h_temp += pop[j]
            e_pop[k] = e_temp
            h_pop[k] = h_temp    
            #loc_e += c_e
        IPR_e = frac_ct**2/sum(e_pop**2)
        IPR_h = frac_ct**2/sum(h_pop**2)
        return [IPR_ex, IPR_e, IPR_h]
