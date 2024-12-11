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

#electronic parameters
class Parameters:
    def __init__(self, Temp, j0, r0j, E_singlet, 
               Epeak, 
               reorE_outer, reorE_inner,
               recom = 0):
        self.Temp, self.Epeak = Temp, Epeak
        self.reorE_outer, self.reorE_inner = reorE_outer, reorE_inner
        self.kT = sp.k/sp.e*self.Temp
        self.j0, self.r0j = j0, r0j
        self.E_singlet = E_singlet
        self.recom = recom    
    def to_dict(self):
        return{
         'Temperature K': self.Temp,
         'Kout s-1ljm ':self.kout,
         'J0_eV':self.j0,
         'R0J_A':self.r0j,
         'E_singlet':self.E_singlet,
         'Epeak':self.Epeak,
         'reorE_inner':self.reorE_inner,
         'reorE_outer':self.reorE_outer         
         }

class site:
    def __init__(self, coordinate, HOMO, LUMO, ID, transition_dipole_Ex, transition_dipole_CT, 
                 recom = 0, Krec_CT = 0, Krec_Ex = 0, V_CT = 0, V_ex = 0):
        self.coordinate = coordinate #XYZ coordiante in space in A
        self.HOMO = HOMO
        self.LUMO = LUMO
        self.Nearest_Neighbour = []
        self.LUMO_coupling = [] 
        self.HOMO_coupling = []
        self.Dipole_coupling = []
        self.transition_dipole_Ex = transition_dipole_Ex
        self.transition_dipole_CT = transition_dipole_CT
        if recom == 0:
            self.Krec_Ex = Krec_Ex            
            self.Krec_CT = Krec_CT
        elif recom == 1:
            self.V_CT = V_CT
            self.V_ex = V_ex
        self.id = ID
        self.recom = recom
            
    def to_dict(self, nonuniform = 0):
        recom = self.recom
        if recom == 0:
            value_ex = self.Krec_Ex
            value_CT = self.Krec_CT
        elif recom == 1:
            value_ex = self.V_ex
            value_CT = self.V_CT
        labels_ex = ['Krec_Ex_s-1', 'V_Ex_eV']
        labels_CT = ['Krec_CT_s-1', 'V_CT_eV']
        return{
            'coordinate': self.coordinate,
            'HOMO':self.HOMO,
            'LUMO':self.LUMO,
            'Nearest_Neighbour':self.Nearest_Neighbour,
            'LUMO_coupling':self.LUMO_coupling,
            'HOMO_coupling':self.HOMO_coupling,
            'Dipole_coupling':self.Dipole_coupling,
            'recom_mech':self.recom,
            labels_ex[recom]:value_ex,
            labels_CT[recom]:value_CT,
            'transition_dipole_Ex':self.transition_dipole_Ex,
            'transition_dipole_CT':self.transition_dipole_CT,
            'id':self.id,
            }

class lattice:
    def __init__(self):
        self.sites=[]
        self.states=[]   
        
    def generate_uniform(self, size:int, HOMO:int, LUMO:int, 
                         dist_sites_A:int, min_dist_near_neighbour,
                         t0_homo, t0_lumo, 
                         d0, r0d, 
                         CT_abs_prop,
                         Krec_Ex = 0, Krec_CT = 0, 
                         V_ex = 0, V_CT = 0,
                         recom = 0):
        counter_id = 0
        self.recom = recom        
           
        site_vec = np.zeros((size**2, 3))
                
        for x,y in product(range(size), repeat = 2):
            #The 1 here is the transition dipole of the exciton
            #To do with how well it absorbs the photons
            if self.recom == 0:
                self.sites.append(site(np.array([x*dist_sites_A, y*dist_sites_A, 0]),
                                       HOMO, LUMO, counter_id, 1, CT_abs_prop, 
                                       Krec_Ex = Krec_Ex, Krec_CT = Krec_CT))
            elif self.recom == 1: 
                self.sites.append(site(np.array([x*dist_sites_A, y*dist_sites_A, 0]),
                                       HOMO, LUMO, counter_id, 1, CT_abs_prop, 
                                       recom = 1, V_CT = V_CT, V_ex = V_ex))
            site_vec[counter_id,:] = np.array([x*dist_sites_A, y*dist_sites_A, 0])
            counter_id = counter_id + 1  
        self.site_vec = site_vec
        for ii,jj in combinations(range(counter_id),2):
            #Find distance from site ii to site jj
            distance_bet_sites = math.dist(self.sites[ii].coordinate, self.sites[jj].coordinate)
            #Only couple sites which are closer than whatever the nearest neighbour distance is 
            if distance_bet_sites < min_dist_near_neighbour:
                self.sites[ii].Nearest_Neighbour.append(self.sites[jj].id)
                self.sites[jj].Nearest_Neighbour.append(self.sites[ii].id)
                if math.isclose(distance_bet_sites-dist_sites_A, 0, abs_tol = 1e-3):
                    self.sites[ii].LUMO_coupling.append(t0_lumo)
                    self.sites[jj].LUMO_coupling.append(t0_lumo)
                    self.sites[jj].HOMO_coupling.append(t0_homo)
                    self.sites[ii].HOMO_coupling.append(t0_homo)
                else:
                    self.sites[ii].LUMO_coupling.append(0)
                    self.sites[jj].LUMO_coupling.append(0)
                    self.sites[jj].HOMO_coupling.append(0)
                    self.sites[ii].HOMO_coupling.append(0)
                self.sites[jj].Dipole_coupling.append(d0*((distance_bet_sites-dist_sites_A)/r0d+1)**-3)
                self.sites[ii].Dipole_coupling.append(d0*((distance_bet_sites-dist_sites_A)/r0d+1)**-3) 

    def Sites_ToDataframe(self, nonuniform = 0):
        return pd.DataFrame.from_records([s.to_dict(nonuniform) for s in self.sites])
    
    #Transdip = transition dipole
    #Set to be one by default for excitons
    #Used as a proxy for strength of generation
    def buildHam(self, params, F, min_dist_near_neighbour, dist_CS_min, disorder_site_ene, random_seed = 0):
        #data is where you store the nergy values which will get put in the Hamiltonian matrix
        #basis is a list of all the possible combinations of electron and hole locations which make up the basis states
        #an element of basis takes the form [Elec_pos.id, Hol_pos.id]
        #Order of basis is to keep the electron id fixed and cycle through all the holes before incrementing electron number by one
        #i.e. (1,1), (1,2), ... (1, n), (2,1)
        #row and col are to keep track of where the elements are in the matrix so that data can be put into the right shape at the end 
        row, col, data, basis = [],[],[],[]
        recom = self.recom
        if random_seed != 0:
            np.random.seed(seed = random_seed)
        #Create empty lists to store the values of various quantities of interest 
        Transdip_vec_Ex, krec_vec_Ex, Transdip_vec_CT, krec_vec_CT, dist_he = [],[],[],[],[]
        E_singlet = params.E_singlet
        Is_Ex, Is_CT = [],[]        
        counter = 0
        for Elec_pos in self.sites:
            for Hol_pos in self.sites: 
                distance_ele_hole = math.dist(Elec_pos.coordinate,Hol_pos.coordinate)
                dist_vec = (Elec_pos.coordinate - Hol_pos.coordinate)
                dist_he.append(distance_ele_hole)
                Coul_j = params.j0/(1+distance_ele_hole/params.r0j)
                #Setting generation and decay rates for CT and exciton states
                #Values set by parameters you input when you build the lattice
                #Don't actually affect the energy calculation but need these for the rates
                if distance_ele_hole == 0:
                    Transdip_vec_Ex.append(Elec_pos.transition_dipole_Ex)
                    Is_Ex.append(1)
                    Is_CT.append(0)
                    #need to do this next 
                    if recom == 0:
                        krec_vec_Ex.append(Elec_pos.Krec_Ex)
                    elif recom == 1:
                        krec_vec_Ex.append(Elec_pos.V_ex)                     
                else:
                    Transdip_vec_Ex.append(0)
                    Is_Ex.append(0)
                    Is_CT.append(1)
                    krec_vec_Ex.append(0)
                if Hol_pos.id in Elec_pos.Nearest_Neighbour:
                    Transdip_vec_CT.append(Elec_pos.transition_dipole_CT)
                    if recom == 0:
                        krec_vec_CT.append(Elec_pos.Krec_CT)
                    #This is making the assumption that the coupling for CT/polaron states which aren't nearest neighbours
                    #to the ground is negligible. 
                    elif recom == 1:
                        krec_vec_CT.append(Elec_pos.V_CT)
                else:
                    Transdip_vec_CT.append(0)
                    krec_vec_CT.append(0)
                    
                #Calculate the values of the diagonal terms which go into the matrix
                #Which diagonal term you are calculating is kept track of my row and col
                #Also build up all the possible combinations of electrons and holes which make up the basis 
                #Keep track of these using basis   
                row.append(counter)
                col.append(counter)
                #Disorder term added on as a random offset in the bandgap energy
                #Note that the sign of the energy contribution will depend on the orientation of the elctron and hole relative to the field
                #One orientation (with hole closer to e.g., LHS than e-) will the high energy and the other (with hole closer to e.g. RHS than e-) will be low energy
                #Should there be disorder in the couplings too?
                if distance_ele_hole == 0:
                    data.append(E_singlet + disorder_site_ene*np.random.standard_normal(1)[0])                     
                else:
                    data.append(Elec_pos.LUMO - Hol_pos.HOMO + disorder_site_ene*np.random.standard_normal(1)[0] - Coul_j + np.dot(dist_vec, F))
                counter += 1
                basis.append([Elec_pos.id, Hol_pos.id])
        
        #off diagonal elements 
        #using combinations so will never have case where ii = jj
        for ii,jj in combinations(range(len(basis)),2):
            #CT state/exciton on basis[ii] and basis[jj] (exciton can only be on one of the two sites)
            if basis[jj][0] in self.sites[basis[ii][0]].Nearest_Neighbour and basis[ii][1] == basis[jj][1]:
                #index() returns the index of the first occurance of the specified value 
                counter_neighbour_ele = self.sites[basis[ii][0]].Nearest_Neighbour.index(basis[jj][0])                        
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
            if basis[ii][0] == basis[ii][1] and basis[jj][0] == basis[jj][1] and basis[jj][0] in self.sites[basis[ii][0]].Nearest_Neighbour:
                counter_neighbour_ele = self.sites[basis[ii][1]].Nearest_Neighbour.index(basis[jj][1])
                row.append(ii)
                col.append(jj)
                data.append(self.sites[basis[ii][0]].Dipole_coupling[counter_neighbour_ele])
                row.append(jj)
                col.append(ii)
                data.append(self.sites[basis[ii][0]].Dipole_coupling[counter_neighbour_ele])
            #CT state/exciton on basis[ii] and basis[jj] (exciton can only be on one of the two sites) 
            if basis[jj][1] in self.sites[basis[ii][1]].Nearest_Neighbour and basis[ii][0]==basis[jj][0]:  
                counter_neighbour_hole = self.sites[basis[ii][1]].Nearest_Neighbour.index(basis[jj][1])                      
                row.append(ii)
                col.append(jj)
                data.append(self.sites[basis[ii][1]].HOMO_coupling[counter_neighbour_hole])
                row.append(jj)
                col.append(ii)
                data.append(self.sites[basis[ii][1]].HOMO_coupling[counter_neighbour_hole])
                        
        #I have made the Hamiltonian unsparse as I worry about numerical errors when using the sparse package         
        #csr_array comes from the scipy.sparse package, builds the sparse array from the data 
        #Do this even if you don't use the sparse method to find eigenvalues as it reduces the memory needed to store the data
        self.Ham = csr_array((data, (row, col)), shape=(len(self.sites)**2,len(self.sites)**2))
        #add more properties to the lattice to be used in rate calculations 
        self.basis = basis
        self.Transdip_vec_CT = Transdip_vec_CT
        self.Transdip_vec_Ex = Transdip_vec_Ex
        self.krec_vec_Ex = krec_vec_Ex
        self.krec_vec_CT = krec_vec_CT         
        self.dist_he = dist_he
        self.Is_Ex = Is_Ex
        self.Is_CT = Is_CT

    #Find the eigenvalues and eigenvectors of the Hamiltonian 
    def states_from_Ham_v2(self, params, max_energy_diff):  
        reorE_inner, reorE_outer = params.reorE_inner, params.reorE_outer, 
        evals,evecs = linalg.eigh(self.Ham.toarray())
        #evals,evecs = eigsh(self.Ham, k=(len(self.sites)**2)-1, which = 'SM')
        evals = np.real(evals)
        evecs = np.real(evecs)
        #remove eigenvectors corresponding to states above the energy cut off
        evecs = evecs[:,evals<evals.min() + max_energy_diff]
        #remove eigenvalues corresponding to states above the energy cut off
        evals = evals[evals<evals.min() + max_energy_diff]
        #Declare a whole ton of empty lists which will be populated by the function 
        list_evecs, dis_st, IPR, Ex_char, Transdip_Ex, Transdip_CT, occupation_prob, Krec_Ex, Krec_CT = [],[],[],[],[],[],[],[],[]
        Is_Ex = np.array(self.Is_Ex)
        basis = np.array(self.basis)
        if params.recom == 1:
            V_eff_ex, V_eff_CT = [],[]
        for i in range(len(evals)):
            list_evecs.append(evecs[:,i])
            #evecs is a list with each element telling you the cn coefficient for the eigenstate in terms of the original set of
            #eignevectors
            #sqaure each cn coefficent - tells you how much of each of the original bais states is in the eigenstate
            occupation_probability = evecs[:,i]**2
            occupation_prob.append(occupation_probability)
            #@ in the middle of a line is a matrix multiplication 
            #occupation_probability is a column vector containing the coefficients of the eignestate 
            #essentailly finding the expectation value of each varaible for the eigenvector under consideration
            #self.dist_he refers tot he separation of the electron and hole in the basis state so what we plot is a weighted sum of these 
            #by the coefficient of that basis state in the eigenstate (i.e., the expectation value)
            dis_st.append(occupation_probability @ self.dist_he)
            #Transition dipole moment used as a proxy for strength of generation
            Transdip_Ex.append(occupation_probability @ self.Transdip_vec_Ex)
            Transdip_CT.append(occupation_probability @ self.Transdip_vec_CT)
            Ex_char.append(occupation_probability @ Is_Ex)
            if params.recom == 0:
                IPR.append(calc_IPR(evecs[:,i], basis, Is_Ex))
                Krec_Ex.append(occupation_probability @ self.krec_vec_Ex)
                Krec_CT.append(occupation_probability @ self.krec_vec_CT)
            #order of vectors is coupling, inner reorgansisation energy, outer reorganisation energy
            #Using generalised MLJ following Taylor and Kassal 2018/ D'Avino et al. J. Phys. Chem. Lett. 2016, 7, 536-540
            elif params.recom == 1:
                IPR.append(calc_IPR(evecs[:,i], basis, Is_Ex))
                Effective_coupling_ex = evecs[:,i]@self.krec_vec_Ex
                V_eff_ex.append(Effective_coupling_ex)
                Effective_inner_reorE_ex = reorE_inner*(1/IPR[i][0])
                Effective_outer_reorE_ex = reorE_outer*(1/IPR[i][0])
                Krec_Ex.append(Recombination.decay_rate(Effective_inner_reorE_ex, params.Epeak, Effective_outer_reorE_ex, 
                                                         evals[i], Effective_coupling_ex))

                Effective_coupling_CT = evecs[:,i]@self.krec_vec_CT
                V_eff_CT.append(Effective_coupling_CT)
                Effective_inner_reorE_CT = reorE_inner*(1/IPR[i][1] + 1/IPR[i][2])
                Effective_outer_reorE_CT = reorE_outer*(1/IPR[i][1] + 1/IPR[i][2])
                Krec_CT.append(Recombination.decay_rate(Effective_inner_reorE_CT, params.Epeak, Effective_outer_reorE_CT, 
                                                         evals[i], Effective_coupling_CT))                                
        states = pd.DataFrame({'En_eV': evals,'dis_eh': dis_st, 
                                'IPR': IPR, 'Ex_char': Ex_char,
                                'Transdip_Ex': Transdip_Ex,'Transdip_CT': Transdip_CT,
                                'Krec_Ex': Krec_Ex,'Krec_CT': Krec_CT,
                                'occupation_probability': occupation_prob})          
        #make it so the states are numbered 1 to n, not 0 to n-1
        states = states.sort_values(by=['En_eV'])
        states.insert(0, "state", np.arange(1, len(evals)+1))
        if params.recom == 1:
            states = states.assign(V_effective_ex = V_eff_ex, V_effective_CT = V_eff_CT)
        if len(self.states) == 0:
            self.states = states   
        else:
            #adding newly calculated values onto the pre-existing states dataframe 
            self.states = pd.concat([self.states, states], ignore_index=True)  
        return np.asarray(list_evecs)
    
    def rates_mat_2D_v7(self, params):
        '''Redfield rates using a secular approximation as in Marcus and Renger 2002 & Quantum biology revisited 2020'''
        Epeak, kT = params.Epeak, params.kT
        reorE_outer, reorE_inner, = params.reorE_outer, params.reorE_inner
        reorE_inner = Redfield.norm_J_Renger(reorE_outer+reorE_inner, Epeak)
        En_eV = np.array(self.states.En_eV)
        num_states = len(En_eV)
        occupation_prob = np.array([self.states.occupation_probability[i] for i in range(num_states)])        
        rates_mat = np.zeros((num_states, num_states), dtype = np.float32)
        inds = np.tril_indices_from(rates_mat, k=-1)
        len_inds = int(((num_states)*(num_states-1))/2) 
        w = np.array([En_eV[inds[0][k]] - En_eV[inds[1][k]] for k in range(len_inds)])
        C = Redfield.C_re_2D_array(w, reorE_inner, Epeak, kT) 
        #Here I pick the kth element from inds[0], which gives the row index, and the kth element from inds[1], which gives the column index
        #This corresponds to a transition between the inds[0][k] and inds[1][k] eigenstates
        gamma_Cw_Re_tot = np.array([occupation_prob[inds[0][k]]@occupation_prob[inds[1][k]]*C[k] for k in range(len_inds)])
        rates_mat[inds] = gamma_Cw_Re_tot
        rates_mat[(inds[1], inds[0])] = rates_mat[inds]*np.exp(-w/params.kT)  
        self.rates_mat = rates_mat

    def solve_steady(self,params):
        Isteady_n = np.zeros(len(self.states))
        #A deep copy is only relevant for a compound objects i.e, one which contains other objects, such as lists
        #In a shallow copy, a new compound object is made and then (to the extent possible) references to the objects in the 
        #original are inserted into it
        #In a deep copy, the new compound object is populated with copies of the object found in the original
        #Advantageous as this lets you make a new object which you can alter without affecting the data stored in the original
        a = deepcopy(self.rates_mat)
        #Originally a[i,j] means the transition from state j to state i
        #Transpose means a[i,j] is now the transition from state i to state j
        a = np.transpose(a)
        self.states['Gen'] = self.states.Transdip_Ex + self.states.Transdip_CT
        
        for i in range(len(self.states)):
            #The sum of a column of the rate matrix gives the total rate of transitions out of the state i (as he's transposed it)
            #Originally a[i,i] is zero as no transitions from a state to itself
            a[i,i] = -sum(a[:,i])
            #this is the total rate of transitions out of i due to both Redfield and other kinetic processes (i.e, recombination, extraction)
            #These diagonal terms will be multiplied by the population of the state i when you do the rate eqn matrix multiplication
            #The other terms in the row will be multiplied by the populations of the other states, which is what you want as this calculates the 
            #rate of transitions from those states into state i
            a[i,i] = a[i,i] - self.states.iloc[i].Krec_Ex - self.states.iloc[i].Krec_CT
            Isteady_n[i] = self.states.iloc[i].Gen
            
        #print(np.linalg.cond(a))
        b = Isteady_n*(-1)
        #linalg.solve computes the “exact” solution, x, of the well-determined, i.e., full rank, linear matrix equation ax = b
        #Here, a are the rates, b is the generation term and x is the SS occupation 
        #So I think z is the SS occupation of the eigenstates of the unperturbed Hamiltonian?
        #Something proportional to this anyway - don't think they've been normalised 
        #I have changed to using linalg.lstsq as this also works for underdetermined matrices, which I sometimes had in the bilayer case
        #I changed back as this one at least gives me warnings that the answers might be wrong...
        #z, res, r, s = linalg.lstsq(a, b)
        z = linalg.solve(a,b)
        self.states['z'] = z
    
#Code needs optimising 
def calc_IPR(eigvec, basis, Is_Ex):
        #These are calculated following method in D'Avino et al. J. Phys. Chem. Lett. 2016, 7, 536-540
        #Although I think there is a mistake in this paper - missing a squared sign???
        #See workbook 3, 14/11/23
        #Very slow currently as there are two for loops
        #Three if you count the for loop this is nested inside in get states
        num_states = len(eigvec)
        num_sites = int(basis[-1][0]+1)
        pop = eigvec**2
        
        frac_FE = sum(pop*Is_Ex)
        frac_CT = sum(pop*(1-Is_Ex))
        IPR_ex = frac_FE**2/sum((pop*Is_Ex)**2)

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
        IPR_e = frac_CT**2/sum(e_pop**2)
        IPR_h = frac_CT**2/sum(h_pop**2)
        return [IPR_ex, IPR_e, IPR_h]
