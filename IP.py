import scipy.interpolate as interpolate
import interp_funcs as interp_funcs
import numpy as np

########################
# RAW ABSORPTION MODEL #
########################
def layer_attenuation(E,layer_thickness,layer_density,element_names,ratios,atomic_weights): # layer_thickness in microns, density in g/cc

    attenuation_data_path = './Attenuation_Data/'
    total_weight = np.sum(np.array(ratios)*np.array(atomic_weights))
    layer_thickness *= 1e-4 # convert from microns to cm

    class element:
        global attenuation_data_path
        global total_weight
        global layer_density
        global mu_en_rho
        def __init__(self, name, ratio, atomic_weight):
            self.name = name
            self.ratio = ratio
            self.atomic_weight = atomic_weight

            attenData          = np.loadtxt(attenuation_data_path+name+'.dat',comments="#")
            photon_MeV, mu_rho = attenData[:,0], attenData[:,1]

            filename = name+'_interp.dat'
            file = open(attenuation_data_path+filename, 'r')
            interp_edge = []
            interp_kind = []
            for line in file:
                p = line.split()
                interp_edge.append(int(p[0]))
                interp_kind.append(p[1])
            file.close()

            if(interp_edge[0]==0):
                self.interp = interp_funcs.no_edge(E,photon_MeV,mu_rho,interp_edge,interp_kind)
            elif(interp_edge[0]==1):
                self.interp = interp_funcs.one_edge(E,photon_MeV,mu_rho,interp_edge,interp_kind)
            elif(interp_edge[0]==2):
                self.interp = interp_funcs.two_edge(E,photon_MeV,mu_rho,interp_edge,interp_kind)
            elif(interp_edge[0]==4):
                self.interp = interp_funcs.four_edge(E,photon_MeV,mu_rho,interp_edge,interp_kind)

            """Plotting Check"""
            if(0):
                import matplotlib.pyplot as plt
                plt.title(self.name)
                plt.plot(photon_MeV*1e3, mu_rho,label='NIST')
                plt.plot(E*1e3,self.interp,label='Interp')
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel('Photon Energy (keV)')
                plt.ylabel(r'$\mu /\rho$')
                plt.legend(loc='best')
                plt.show()

            self.density = ratio * (atomic_weight/total_weight) * layer_density
            self.mass_thickness = layer_thickness * self.density
            self.attenuation_factor = np.exp(-(self.interp*self.mass_thickness))

    elements = list()
    for element_name,ratio,atomic_weight in zip(element_names,ratios,atomic_weights):
        elements.append(element(element_name,ratio,atomic_weight))

    #plotting check
    if(0):
        import matplotlib.pyplot as plt
        factor = 1.
        for element in elements:
            factor *= element.interp
            plt.plot(E*1e3,element.interp,label=element.name)
        # plt.plot(E*1e3,factor,label='total')
        plt.legend(loc='best')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Photon Energy (keV)')
        plt.ylabel(r'$\mu /\rho$')
        plt.show()

    attenuation_factor = 1.
    for element in elements:
        attenuation_factor *= element.attenuation_factor

    return attenuation_factor

#########
#FITTING#
#########

# Fitting in needed to match the raw absoprtion model to the experimental data

def fitting_function(m,E,b):
    return m*E + b

def scaling(E,type):
    #Type decides between Meadowcroft (True) and Maddox (False)
    E *= 1e6

    booles_E_I = (E<(5.98880e3))
    inds_E_I   = np.where(booles_E_I)
    E_I        = E[inds_E_I]

    booles_E_II = np.logical_and(E>=5.98880e3,E<13.4737e3)
    inds_E_II   = np.where(booles_E_II)
    E_II        = E[inds_E_II]

    booles_E_III = np.logical_and(E>=13.4737e3,E<37.4406e3)
    inds_E_III   = np.where(booles_E_III)
    E_III        = E[inds_E_III]

    booles_E_IV = (E>=(37.4406e3))
    inds_E_IV   = np.where(booles_E_IV)
    E_IV        = E[inds_E_IV]


    if(type):
        #Meadowcroft
        MR_Region_I_m = 1.07*1e-3
        MR_Region_I_b = -0.51

        MR_Region_II_m = 0.61*1e-3
        MR_Region_II_b = 1.92

        MR_Region_III_m = 0.67*1e-3
        MR_Region_III_b = 0.75

        MR_Region_IV_m = 0.67*1e-3
        MR_Region_IV_b = -9.10

        # plt.title('Meadowcroft Sensitivity')

        Sensitivity_MR_scaling = np.concatenate((fitting_function(MR_Region_I_m, E_I, MR_Region_I_b),fitting_function(MR_Region_II_m, E_II, MR_Region_II_b)))
        Sensitivity_MR_scaling = np.concatenate((Sensitivity_MR_scaling,fitting_function(MR_Region_III_m, E_III, MR_Region_III_b)))
        Sensitivity_MR_scaling = np.concatenate((Sensitivity_MR_scaling,fitting_function(MR_Region_IV_m, E_IV, MR_Region_IV_b)))

    else:
        #Maddox
        MR_Region_II_m = 7.098e-4
        MR_Region_II_b = -3.789

        MR_Region_III_m = 4.446e-4
        MR_Region_III_b = 0.3745

        MR_Region_IV_m = 5.658e-4
        MR_Region_IV_b = -12.05

        # plt.title('Maddox Sensitivity')

        Sensitivity_MR_scaling = np.concatenate((fitting_function(MR_Region_II_m, E_II, MR_Region_II_b),fitting_function(MR_Region_III_m, E_III, MR_Region_III_b)))
        Sensitivity_MR_scaling = np.concatenate((Sensitivity_MR_scaling,fitting_function(MR_Region_IV_m, E_IV, MR_Region_IV_b)))

    mPSL_scaling = Sensitivity_MR_scaling
    PSL_scaling  = mPSL_scaling*1e-3
    E *= 1e-6

    return PSL_scaling

def sensitivity(E,phosphor_width,phosphor_density,phosphor_elements,phosphor_ratios,phosphor_atomic_weights,type=True):

    attenuation = layer_attenuation(E,phosphor_width,phosphor_density,phosphor_elements,phosphor_ratios,phosphor_atomic_weights)
    absorption = 1. - attenuation

    scaling_factor = scaling(E,type)

    sens = absorption * scaling_factor

    """
    testing
    """
    if(0):
        import matplotlib.pyplot as plt
        plt.plot(1e3*E,absorption,label='absorption')
        plt.plot(1e3*E,attenuation,label='attenuation')
        plt.xlabel('Energy (keV)')
        plt.legend(loc='best')
        plt.show()

        # this plot will seem x2 as large as it should be, this is down to the scaling that is used in the next step
        plt.plot(1e3*E,1e3*E*absorption,label='E*abs')
        # np.savetxt('E*abs.txt',1e3*E*absorption,fmt='%.5e')
        plt.xlabel('Energy (keV)')
        plt.legend(loc='best')
        plt.show()

        plt.plot(1e3*E,1e3*absorption*scaling(E,type),label='abs*scaling')
        # np.savetxt('abs*scaling.txt',1e3*E*sens,fmt='%.5e')
        plt.xlabel('Energy (keV)')
        plt.legend(loc='best')
        plt.grid(True)
        plt.show()

        plt.plot(1e3*E,sens,label='sens')
        # np.savetxt('abs*scaling*1_E.txt',sens,fmt='%.5e')
        plt.xlabel('Energy (keV)')
        plt.legend(loc='best')
        plt.grid(True)
        plt.show()

    return sens

#######################
# IP LAYER DEFINITION #
#######################

class layer:
    global atomic_weights_dict
    atomic_weights_dict = {'Bromine':79.9,'Fluorine':19.0,'Barium':137.3,'Iodine':126.9,'Mn':54.9,'Zn':65.4,'Fe':55.8,'C':12.0,'H':1.0,'O':16.0,'Mylar':1.,'Polystyrene':1.}
    def __init__(self, E, name, width, density, elements, ratios, sens_type=True):
        self.E = E # Energy range
        self.name = name # Layer name as string
        self.width = width # in microns
        self.density = density # in g/cc
        self.elements = elements # Element names ans strings
        self.ratios = ratios # Array of element ratios
        self.sens_type = sens_type # Bool for sensitivity model; True for Meadowcroft, False for Maddox
        self.atomic_weights = [atomic_weights_dict[element] for element in self.elements]
        self.attenuation = layer_attenuation(E,self.width,self.density,self.elements,self.ratios,self.atomic_weights)
        if(self.name=='Phosphor'):
            self.sensitivity = sensitivity(E,self.width,self.density,self.elements,self.ratios,self.atomic_weights,self.sens_type)
