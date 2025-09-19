from tc_python import *
import numpy as np
import math
import time
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

class alloy_sys:
    """
    
    """
    def __init__(self,dependent='Al',composition={"Si":0.1,"Mg":0.0035},param="T",param_range=range(300,1300,100),custom=False,database='TCAL9',unit='X',suspended=None,default=True,included=None):
        self.dependent = dependent
        self.composition = composition
        self.param=param
        self.param_range = param_range
        self.unit = unit
        self.custom = custom
        self.database = database
        self.suspended = suspended
        self.default = default
        self.included = included
    def add_phases(self,phaselists):
        self.phases=phaselists
    def add_nominal_compositions(self,nom_comp):
        self.nominal_comp=nom_comp
    def add_phase_dict(self,phase_dict):
        self.phase_dict = phase_dict
    def add_elem_dict(self,elemcomps):
        self.elem_dict = elemcomps
    def do_perform_single_axis_split(self,T=298,threshold=0.0001):
        """
        # This function is intended to be a general script to perform single acix calculations using a local ThermoCalc (TC-Python) installation
        # It is currently configured to accept temperature or a non-balancing element as the "axis"
        # Inputs:
        # "alloy_sys" class that contains dependent balancing element (defualt dependent="Al"), composition of alloying elements (default composition={"Si": 0.1,"Mg": 0.0035}), 
        # element unit ("W" for compositions in weight fraction, "X" for (defualt unit="W")), parameter varied as either temperature or alloying elemenet (default param="T"),
        # range parameter is varied on (default param_range=range(300,1200,100)), a boolean statement for if a custom database is in use (default custom=False),
        # and database used for calculations (as TC database or path to user database, defualt="TCAL9")
        """
        alloy=self
        elements = [alloy.dependent] + list(alloy.composition.keys())
        phaselists = []
        nom_comps = []
        if alloy.param !='T' and alloy.param not in elements:
            raise ValueError("Invalid parameter [%s] selected, Use Temperature [T] or a non-balancing element"%param+str(list(composition.keys()))) 
        if (alloy.default == False) & (alloy.included is None):
            raise ValueError("No phases defined - Default phases disabled but no phases selected")
        with TCPython(logging_policy=LoggingPolicy.NONE)as start:
            if alloy.custom:
                if alloy.default:
                    gibbs = start.select_user_database_and_elements(alloy.database,elements).get_system();
                else:
                    gibbs = (start.select_user_database_and_elements(alloy.database,elements)
                            .without_default_phases())
                    for ph in alloy.included:
                        gibbs = gibbs.select_phase(ph)
                    gibbs = gibbs.get_system()
            else:
                if alloy.default:
                    gibbs = start.select_database_and_elements(alloy.database,elements).get_system();
                else:
                    gibbs = (start.select_database_and_elements(alloy.database,elements)
                            .without_default_phases())
                    for ph in alloy.included:
                        gibbs = gibbs.select_phase(ph)
                    gibbs = gibbs.get_system()
            eq_calculation = gibbs.with_single_equilibrium_calculation();            
            for r in alloy.param_range:
                temp_comp = {}
                if alloy.param in elements:
                    alloy.composition[alloy.param] = r
                if alloy.param == 'T':
                    T=r
                eq_calculation.set_condition("T",T)
                tot = 0
                for e in list(alloy.composition.keys()):
                    eq_calculation.set_condition("%s(%s)"%(alloy.unit,e),alloy.composition[e])
                    temp_comp[e] = alloy.composition[e]
                    tot +=alloy.composition[e]
                temp_comp[alloy.dependent] = 1-tot
                if alloy.suspended is not None:
                    for phase in alloy.suspended:
                        eq_calculation = eq_calculation.set_phase_to_suspended(phase)                    
                calc_result = eq_calculation.calculate()
                stable_phases = calc_result.get_stable_phases()
                tempdict = {}
                for phase in stable_phases:
                    tar = {}
                    tar['frac'] = calc_result.get_value_of('NP(' + phase + ')')
                    for e in elements:
                        tar[e] = calc_result.get_value_of('X(%s,%s)'%(phase,e))
                    tempdict[phase] = tar
                phaselists.append(tempdict)
                nom_comps.append(temp_comp)
            phases = []
            for p in phaselists:
                phases.extend(list(p.keys()))
            keyphases = np.unique(phases)
            num_phases = len(keyphases)
            totelems = len(elements)
            phase_dict = {}
            for k in keyphases:
                phase_dict[k]={}
                phase_dict[k]['frac'] = np.zeros(len(phaselists))
                for e in elements:
                    phase_dict[k][e] = np.zeros(len(phaselists))
                for p in range(len(phaselists)):
                    if k in phaselists[p]:
                        phase_dict[k]['frac'][p] = phaselists[p][k]['frac']
                        for e in elements:
                            phase_dict[k][e][p] = phaselists[p][k][e]
            elemcomps = {}
            comp_trend = {}
            for e in elements:
                comp_trend[e] = np.zeros(len(phaselists))
                for n in range(len(nom_comps)):
                    comp_trend[e][n] = nom_comps[n][e]
            for e in elements:
                elemcomps[e]={}
                for k in keyphases:
                    elemcomps[e][k]=phase_dict[k][e]*phase_dict[k]['frac']/comp_trend[e]
        alloy.add_phases(phaselists)
        alloy.add_nominal_compositions(comp_trend)
        alloy.add_phase_dict(phase_dict)
        alloy.add_elem_dict(elemcomps)
    def phase_distribution(self,save=False,img_path='PhraseDistribution.png',phase_list=[None],excluded=True,
                           color_list=['black','red','magenta','blue','darkorange',
                                       'darkturquoise','cyan','blueviolet']):
        # old colorlist: color_list=['black','darkgreen','darkorange','cornflowerblue','slategrey','magenta','lime','mediumturquoise','blue']
        alloy=self
        fig,ax = plt.subplots(figsize=(6,5))
        keyphases = list(alloy.phase_dict.keys())
        k = 0
        threshold = 0.0001
        lines='solid'
        # xlab = "Temperature (K)"
        # xparams = list(alloy.param_range)
        xlab = "Temperature (C)"
        xparams = [i-273 for i in list(alloy.param_range)]
        if alloy.param in alloy.composition.keys():
            xlab = alloy.param + " ("+alloy.unit+". % )"
            xparams = alloy.param_range*100
        if excluded:
            for kp in keyphases:
                if any(alloy.phase_dict[kp]['frac']>threshold)&(kp not in phase_list):
                    if k>=len(color_list):
                        k=0
                        lines='dashed'
                    ax.plot(xparams,alloy.phase_dict[kp]['frac']*100,c=color_list[k],label=kp,linestyle=lines,linewidth=5)
                k+=1
        else:
            for kp in keyphases:
                if any(alloy.phase_dict[kp]['frac']>threshold)&(kp in phase_list):
                    if k>=len(color_list):
                        k=0
                        lines='dashed'
                    ax.plot(xparams,alloy.phase_dict[kp]['frac']*100,c=color_list[k],label=kp,linestyle=lines,linewidth=5)
                k+=1
        ax.legend(loc='upper left',bbox_to_anchor=(1,0.85))
        ax.set_xlim((xparams[0],xparams[-1]))
        ax.set_ylim((1e-2,150))#10e2))
        ax.set_yscale('log')
        ax.set_xlabel(xlab,fontweight='bold');
        ax.set_ylabel("Amount of phase (at. %)",fontweight='bold')
        ax.set_title("Phase Distribution",fontweight='bold',fontsize=15)
        if save:
            plt.savefig(img_path,bbox_inches='tight')
        plt.show()
    def composition_distribution(self,save=False,img_path='CompositionDistribution.png',threshold=0.00001,phase_list=[None],excluded=True,
                                 color_list=['black','red','magenta','blue','darkorange',
                                       'darkturquoise','cyan','blueviolet']):
        # old colorlist: color_list=['black','darkgreen','darkorange','cornflowerblue','slategrey','magenta','lime','mediumturquoise','blue']
        alloy=self
        keyphases = list(alloy.phase_dict.keys())
        num_phases = len(keyphases)
        if excluded:
            if phase_list[0] is not None:
                num_phases = len(keyphases) - len(phase_list)
                keyphases = list(set(keyphases)^set(phase_list))
        else:
            keyphases = phase_list
            num_phases = len(phase_list)
        elements = [alloy.dependent] + list(alloy.composition.keys())
        fig,axes = plt.subplots(math.ceil(num_phases/3),3,figsize=(12,3.25*math.ceil(num_phases/3)))
        k = 0
        col = 0
        # xlab = "Temperature (K)"
        # xparams = list(alloy.param_range)
        xlab = "Temperature (C)"
        xparams = [i-273 for i in list(alloy.param_range)]
        elem_fonts = {}
        cid = 0
        lines='solid'
        for i in range(len(elements)):
            if cid >= len(color_list):
                cid = 0
                lines = 'dashed'
            elem_fonts[elements[i]] = [color_list[cid],lines]
            cid+=1
        if alloy.param in alloy.composition.keys():
            xlab = alloy.param + " ("+alloy.unit+". % )"
            xparams = alloy.param_range*100
        if num_phases > 3:
            for i in range(math.ceil(num_phases/3)):
                for j in range(3):
                    if k<num_phases:
                        # if k>=len(elements):
                        #     cid=0
                        #     lines='dashed'
                        kp = keyphases[k]
                        if any(alloy.phase_dict[kp]['frac']>threshold):#&(kp not in excluded):
                            for e in elements:
                                if (any(alloy.phase_dict[kp][e] > 0.00001)):#&(kp not in excluded):
                                    # axes[i][j].plot(xparams,100*alloy.phase_dict[kp][e],c=color_list[elements.index(e)],linestyle=lines)
                                    axes[i][j].plot(xparams,100*alloy.phase_dict[kp][e],c=elem_fonts[e][0],linestyle=elem_fonts[e][1],linewidth=5)
                        axes[i][j].set_title(kp,fontweight='bold',fontsize=15)
                        axes[i][j].set_yscale('log')
                        # axes[i][j].locator_params(axis='y', numticks=8)
                        axes[i][j].set_xlim((xparams[0],xparams[-1]))
                        axes[i][j].set_ylim((1e-2,150))#10e2))
                        # axes[i][j].set_yticks([1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2][1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2])
                        k+=1
                        cid+=1
                        if i == (math.ceil(num_phases/3)-1):
                            axes[i][j].set_xlabel(xlab,fontweight='bold')
                        if j==0:
                            axes[i][j].set_ylabel("Composition of phase (at. %)",fontweight='bold')
                    else:
                        axes[i][j].set_visible(False)
                        axes[i-1][j].set_xlabel(xlab,fontweight='bold')
            lines=[Line2D([0],[0],color=color_list[elements.index(e)]) for e in elements]
            axes[0][0].legend(lines,elements,loc='upper right',bbox_to_anchor=(3.75,0.75))
            fig.subplots_adjust(hspace=0.25)
            # fig.tight_layout()
        else:
            for i in range(3):
                if k<num_phases:
                    if k>=len(elements):
                        cid=0
                        lines='dashed'
                    kp = keyphases[k]
                    if any(alloy.phase_dict[kp]['frac']>threshold):
                        for e in elements:
                            if any(alloy.phase_dict[kp][e] > 0.00001)&(kp not in phase_list):
                                axes[i].plot(xparams,100*alloy.phase_dict[kp][e],c=color_list[elements.index(e)],linestyle=lines,linewidth=5)
                    axes[i].set_title(kp,fontweight='bold',fontsize=15)
                    axes[i].set_yscale('log')
                    axes[i].set_xlabel(xlab,fontweight='bold')
                    # axes[i][j].locator_params(axis='y', numticks=8)
                    axes[i].set_xlim((xparams[0],xparams[-1]))
                    axes[i].set_ylim((1e-2,150))#5e2))
                    # axes[i][j].set_yticks([1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2][1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2])
                    k+=1
                    cid+=1
                else:
                    axes[i].set_visible(False)
            lines=[Line2D([0],[0],color=color_list[elements.index(e)]) for e in elements]
            axes[0].legend(lines,elements,loc='upper right',bbox_to_anchor=(3.75,0.75))
            axes[0].set_ylabel("Composition of phase (at. %)",fontweight='bold')
            # plt.subplots_adjust(top=0.75)
        # plt.suptitle("Composition of Phases",fontweight='bold',fontsize=20)
        if save:
            plt.savefig(img_path,bbox_inches='tight')
        plt.show()
    def element_distribution(self,save=False,img_path='ElementDistribution.png',threshold=0.00001,phase_list=[None],excluded=True,
                             color_list=['black','red','magenta','blue','darkorange',
                                       'darkturquoise','cyan','blueviolet']):
        # old colorlist: color_list=['black','darkgreen','darkorange','cornflowerblue','slategrey','magenta','lime','mediumturquoise','blue']
        alloy=self
        keyphases = list(alloy.phase_dict.keys())
        num_phases = len(keyphases)
        if excluded:
            if phase_list[0] is not None:
                num_phases = len(keyphases) - len(phase_list)
                keyphases = list(set(keyphases)^set(phase_list))
        else:
            keyphases = phase_list
            num_phases = len(keyphases)
        elements = [alloy.dependent] + list(alloy.composition.keys())
        num_elems = len(elements)
        k = 0
        col = 0
        lines=[Line2D([0],[0],color=color_list[list(keyphases).index(k)]) for k in keyphases]
        # xlab = "Temperature (K)"
        # xparams = list(alloy.param_range)
        xlab = "Temperature (C)"
        xparams = [i-273 for i in list(alloy.param_range)]
        if alloy.param in alloy.composition.keys():
            xlab = alloy.param + " ("+alloy.unit+". % )"
            xparams = alloy.param_range*100
        if num_elems>3:
            fig,axes = plt.subplots(math.ceil(num_elems/3),3,figsize=(12,3.25*math.ceil(num_elems/3)))
            fig.subplots_adjust(hspace=0.25)
            for i in range(math.ceil(num_elems/3)):
                for j in range(3):
                    if k<num_elems:
                        e = elements[k]
                        for kp in keyphases:
                            if (any(alloy.elem_dict[e][kp]>threshold)):#&(kp not in excluded):
                                axes[i][j].plot(xparams,100*alloy.elem_dict[e][kp],c=color_list[list(keyphases).index(kp)],linewidth=5)
                        axes[i][j].set_title(e,fontweight='bold',fontsize=15)
                        axes[i][j].set_yscale('log')
                        axes[i][j].set_xlim((xparams[0],xparams[-1]))
                        axes[i][j].set_ylim((1e-2,150))#10e2))
                        k+=1
                        if j == 0:
                            axes[i][j].set_ylabel("Amount in phase (%)",fontweight='bold')
                        if i == math.ceil(num_elems/3) - 1:
                            axes[i][j].set_xlabel(xlab,fontweight='bold')
                    else:
                        axes[i][j].set_visible(False)
                        axes[i-1][j].set_xlabel(xlab,fontweight='bold')
                axes[0][0].legend(lines,keyphases,loc='upper right',bbox_to_anchor=(4.1,0.75))
                # fig.tight_layout()
                # fig.subplots_adjust(wspace=1)
        else:
            fig,axes = plt.subplots(1,3,figsize=(12,3.25))
            for i in range(num_elems):
                if k<num_elems:
                    e = elements[k]
                    for kp in keyphases:
                        if any(alloy.elem_dict[e][kp]>threshold):#&(kp not in excluded):
                            axes[i].plot(xparams,100*alloy.elem_dict[e][kp],c=color_list[list(keyphases).index(kp)],linewidth=5)
                    axes[i].set_title(e,fontweight='bold',fontsize=15)
                    axes[i].set_yscale('log')
                    axes[i].set_xlim((xparams[0],xparams[-1]))
                    axes[i].set_ylim((1e-2,150))#5e2))
                    axes[i].set_xlabel(xlab,fontweight='bold')
                    k+=1
                else:
                    axes[i].set_visible(False)
            axes[0].set_ylabel("Amount in phase (%)",fontweight='bold')
            axes[0].legend(lines,keyphases,loc='upper right',bbox_to_anchor=(4.1,0.75))
            plt.subplots_adjust(top=0.75)
        # plt.suptitle("Element Distribution in Phases",fontweight='bold',fontsize=20)
        if save:
            plt.savefig(img_path,bbox_inches='tight')
        plt.show()