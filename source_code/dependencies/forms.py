# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 15:06:26 2019

@author: ani228
"""

from dolfin import *
import sys
import numpy as np

class Forms(object):

    def __init__(self, params):    # amir: sth like constructor

        self.parameters = self.default_parameters()
        self.parameters.update(params)

    def default_parameters(self):
        return {#"bff"  : 29.0,
			#"bfx"  : 13.3,
			#"bxx"  : 26.6,
			"Kappa": 1e5,
			"incompressible" : True,
			};


    def Fmat(self):

        u = self.parameters["displacement_variable"]
        d = u.ufl_domain().geometric_dimension()
        I = Identity(d)
        F = I + grad(u)
        return F

    def Fe(self):
        Fg = self.parameters["growth_tensor"]
        F = self.Fmat()
        if (Fg is None):
            Fe = F
        else:
            print "calculating Fe"
            Fe = as_tensor(F[i,j]*inv(Fg)[j,k], (i,k))
        return Fe

    def Emat(self):

        u = self.parameters["displacement_variable"]
        d = u.ufl_domain().geometric_dimension()
        I = Identity(d)
        #F = self.Fmat()
    	F = self.Fe()
        #return 0.5*(F.T*F-I)
    	return 0.5*(as_tensor(F[k,i]*F[k,j] - I[i,j], (i,j)))


    def Cmat(self):

        u = self.parameters["displacement_variable"]
        d = u.ufl_domain().geometric_dimension()
        #F = self.Fmat()
        F = self.Fe()
        #return F.T*F
        return as_tensor(F[k,i]*F[k,j],(i,j))

    def J(self):
        #F = self.Fmat()
        F = self.Fe()
        return det(F)


    def LVcavityvol(self):

        u = self.parameters["displacement_variable"]
        N = self.parameters["facet_normal"]
        mesh = self.parameters["mesh"]
        X = SpatialCoordinate(mesh)
        ds = dolfin.ds(subdomain_data = self.parameters["facetboundaries"])

        F = self.Fmat()

        vol_form = -Constant(1.0/3.0) * inner(det(F)*dot(inv(F).T, N), X + u)*ds(self.parameters["LVendoid"])

        return assemble(vol_form, form_compiler_parameters={"representation":"uflacs"})

    def RVcavityvol(self):

        u = self.parameters["displacement_variable"]
        N = self.parameters["facet_normal"]
        mesh = self.parameters["mesh"]
        X = SpatialCoordinate(mesh)
        ds = dolfin.ds(subdomain_data = self.parameters["facetboundaries"])

        F = self.Fmat()

        vol_form = -Constant(1.0/3.0) * inner(det(F)*dot(inv(F).T, N), X + u)*ds(self.parameters["RVendoid"])

        return assemble(vol_form, form_compiler_parameters={"representation":"uflacs"})


    def LVcavitypressure(self):

        W = self.parameters["mixedfunctionspace"]
        w = self.parameters["mixedfunction"]
        mesh = self.parameters["mesh"]

        comm = W.mesh().mpi_comm()
        dofmap =  W.sub(self.parameters["LVendo_comp"]).dofmap()
        val_dof = dofmap.cell_dofs(0)[0]

	    # the owner of the dof broadcasts the value
        own_range = dofmap.ownership_range()

        try:
            val_local = w.vector()[val_dof][0]
        except IndexError:
                val_local = 0.0


        pressure = MPI.sum(comm, val_local)

        return pressure



    def RVcavitypressure(self):

        W = self.parameters["mixedfunctionspace"]
        w = self.parameters["mixedfunction"]
        mesh = self.parameters["mesh"]

        comm = W.mesh().mpi_comm()
        dofmap =  W.sub(self.parameters["RVendo_comp"]).dofmap()
        val_dof = dofmap.cell_dofs(0)[0]

	    # the owner of the dof broadcasts the value
        own_range = dofmap.ownership_range()

        try:
            val_local = w.vector()[val_dof][0]
        except IndexError:
            val_local = 0.0


        pressure = MPI.sum(comm, val_local)

        return pressure

    def TempActiveStress(self,time):
        print "inside tempactivestress"
        f0 = self.parameters["fiber"]
        cbforce = Expression('A*(B+sin((B/C)*time + D))', A=15000., B=1., C=16., D=80.2, time = time, degree=0)

        Pactive = cbforce * as_tensor(f0[i]*f0[j], (i,j))
        return Pactive, cbforce

    def PassiveMatSEF(self,hsl):
    #def PassiveMatSEF(self):
        Ea = self.Emat()
        f0 = self.parameters["fiber"]
        s0 = self.parameters["sheet"]
        n0 = self.parameters["sheet-normal"]
        C2 = self.parameters["c2"]
        C3 = self.parameters["c3"]
        bff = self.parameters["bf"][0]
        bfx = self.parameters["bt"][0]
        bxx = self.parameters["bfs"][0]
        Kappa = self.parameters["Kappa"]
        isincomp = self.parameters["incompressible"]
        Cmat = self.Cmat()
        phi_m = self.parameters["phi_m"][0]
        phi_g = self.parameters["phi_g"][0]
        hsl0 = self.parameters["hsl0"]

        if(isincomp):
            p = self.parameters["pressure_variable"]

        C = self.parameters["c"]

        Eff = inner(f0, Ea*f0)
        Ess = inner(s0, Ea*s0)
        Enn = inner(n0, Ea*n0)
        Efs = inner(f0, Ea*s0)
        Efn = inner(f0, Ea*n0)
        Ens = inner(n0, Ea*s0)

        alpha = sqrt(2.0 * Eff + 1.0)
        myofiber_stretch = hsl/hsl0
        print "Myofiber stretch:"
        print project(myofiber_stretch,FunctionSpace(self.parameters["mesh"],"DG",0)).vector().array()
        #alpha = sqrt(dot(f0, Cmat*f0))


    		#QQ = bff*pow(Eff,2.0) + bfx*(pow(Ess,2.0)+ pow(Enn,2.0)+ 2.0*pow(Ens,2.0)) + bxx*(2.0*pow(Efs,2.0) + 2.0*pow(Efn,2.0))

        #QQ_m = conditional(alpha > 1.0, C3*(alpha - 1.0)**2.0, 0.0)
        QQ_m = conditional(myofiber_stretch > 1.0, C3*(myofiber_stretch - 1.0)**2.0, 0.0)
        #QQ_m = C3*(alpha - 1.0)**2.0

        QQ_c = bff*Eff**2.0 + bfx*(Ess**2.0 + Enn**2.0 + 2.0*Ens**2.0) + bxx*(2.0*Efs**2.0 + 2.0*Efn**2.0)
        #QQ_i = (C/2)*Eff**2 + bfx*(Ess**2.0 + Enn**2.0 + 2.0*Ens**2.0) + bxx*(2.0*Efs**2.0 + 2.0*Efn**2.0)



        Wp_m = C2*(exp(QQ_m) -  1.0)

        Wp_m_weighted = phi_m*Wp_m

        if(isincomp):
            Wp_c = C/2.0*(exp(QQ_c) -  1.0) - p*(self.J() - 1.0)
        else:
            Wp_c = C/2.0*(exp(QQ_c) -  1.0) + Kappa/2.0*(self.J() - 1.0)**2.0

        Wp_c_weighted = phi_g*Wp_c

        Wp = Wp_m + Wp_c
        #Wp = Wp_m_weighted + Wp_c_weighted
        #Wp = Wp_c
        return Wp


    def LVV0constrainedE(self):


        mesh = self.parameters["mesh"]
        u = self.parameters["displacement_variable"]
        ds = dolfin.ds(subdomain_data = self.parameters["facetboundaries"])
        dsendo = ds(self.parameters["LVendoid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"])
        pendo = self.parameters["lv_volconst_variable"]
        V0= self.parameters["lv_constrained_vol"]

        X = SpatialCoordinate(mesh)
        x = u + X

        F = self.Fmat()
        N = self.parameters["facet_normal"]
        n = cofac(F)*N

        area = assemble(Constant(1.0) * dsendo, form_compiler_parameters={"representation":"uflacs"})
        V_u = - Constant(1.0/3.0) * inner(x, n)
        Wvol = (Constant(1.0/area) * pendo  * V0 * dsendo) - (pendo * V_u *dsendo)

        return Wvol


    def RVV0constrainedE(self):


        mesh = self.parameters["mesh"]
        self.parameters["displacement_variable"]
        ds = dolfin.ds(subdomain_data = self.parameters["facetboundaries"])
        dsendo = ds(self.parameters["RVendoid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"])
        pendo = self.parameters["rv_volconst_variable"]
        V0= self.parameters["rv_constrained_vol"]

        X = SpatialCoordinate(mesh)
        x = u + X

        F = self.Fmat()
        N = self.parameters["facet_normal"]
        n = cofac(F)*N

        area = assemble(Constant(1.0) * dsendo, form_compiler_parameters={"representation":"uflacs"})
        V_u = - Constant(1.0/3.0) * inner(x, n)
        Wvol = (Constant(1.0/area) * pendo  * V0 * dsendo) - (pendo * V_u *dsendo)

        return Wvol


# this returns only myscle stress (no collagen contribution)
    def stress(self,hsl):
    #def stress(self):
        #set_point = 5000 # kurtis trying something
        mesh = self.parameters["mesh"]

        e1 = Constant((1.0, 0.0, 0.0))
        e2 = Constant((0.0, 1.0, 0.0))
        e3 = Constant((0.0, 0.0, 1.0))

        C = self.parameters["c"]
        bff = self.parameters["bf"][0]
        bfx = self.parameters["bt"][0]
        bxx = self.parameters["bfs"][0]

        f0 = self.parameters["fiber"]
        s0 = self.parameters["sheet"]
        n0 = self.parameters["sheet-normal"]

        C2 = self.parameters["c2"]
        C3 = self.parameters["c3"]
        phi_m = self.parameters["phi_m"][0]
        phi_g = self.parameters["phi_g"][0]
        isincomp = self.parameters["incompressible"]
        hsl0 = self.parameters["hsl0"]

        if(isincomp):
            p = self.parameters["pressure_variable"]
        u = self.parameters["displacement_variable"]
        d = u.ufl_domain().geometric_dimension()
        I = Identity(d)
        #F = self.Fmat()
        F = self.Fe()
        #F=Fe
        J = self.J()
        Ea = self.Emat()

        # Make it a variable to try to differentiate wrt?
        Ea = variable(Ea)

        Eff = inner(f0, Ea*f0)
        Ess = inner(s0, Ea*s0)
        Enn = inner(n0, Ea*n0)
        Efs = inner(f0, Ea*s0)
        Efn = inner(f0, Ea*n0)
        Ens = inner(n0, Ea*s0)

        alpha = sqrt(2.0 * Eff + 1.0)
        myofiber_stretch = hsl/hsl0


        #Q = C3*conditional(alpha>1.0,alpha - 1.0,0.0)**2.0
        ### KURTIS START HERE, SFF NEEDS TO BE ZERO FOR MYOFIBER_STRETCH < 1 ****
        Q = C3*conditional(myofiber_stretch > 1.0, myofiber_stretch - 1.0,0.0)**2.0
        Sff = 2.0 * C2 * C3 * (1.0 - conditional(myofiber_stretch > 1.0,1.0/myofiber_stretch,1.0)) * exp(Q)
        #Sff = 2.0 * C2 * C3 * (1.0 - conditional(alpha > 1.0,1.0/alpha,1.0)) * exp(Q)
        Sff_weighted = Sff*phi_m

        # Calculate Guccione passive stress?
        Qbulk = bff*Eff**2.0 + bfx*(Ess**2.0 + Enn**2.0 + 2.0*Ens**2.0) + bxx*(2.0*Efs**2.0 + 2.0*Efn**2.0)
        #test_E = Function(Quad)
        #test_E.vector()[:] = 1.01
        #Qbulk_test = bff*test_E**2.0
        #Wp_c_test = C/2.0*(exp(Qbulk_test) - 1.0) - p*(self.J() - 1.0)
        #sbulk_test = diff(Wp_c_test,Ea)
        #Pg_test = F*sbulk_test
        #test_passive_stress = inner(f0,Pg_test*f0)


        Wp_c = C/2.0*(exp(Qbulk) -  1.0) - p*(self.J() - 1.0)
        Wp_c_weighted = Wp_c*phi_g

        #sbulk differentiated wrt Ea, in global coordinates
        #thus sbulk is PK2 in global?
        sbulk = diff(Wp_c_weighted,Ea)
        #eff_final = test_E
        #test_value = project(test_passive_stress-set_point, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        #test_v = interpolate(test_value,Quad)
        #test_v_array = test_v.vector().get_local()[:]
        #test_max = np.amax(test_v_array)
        #print "test_max = " + str(test_max)
        """if mm > 10:

            while np.abs(test_max) > 10:
                print "in while loop"
                test_E -= test_passive_stress/inner(f0,diff(Pg_test,Ea)*f0)
                #test_E_fcn= project(test_E,FunctionSpace(mesh,"DG",1),form_compiler_parameters={"representation":"uflacs"})
                #print test_E_fcn.vector().get_local()[:]
                Qbulk_test = bff*test_E**2.0
                Wp_c_test = C/2.0*(exp(Qbulk_test) - 1.0) - p*(self.J() - 1.0)
                sbulk_test = diff(Wp_c_test,Ea)
                Pg_test = F*sbulk_test
                test_passive_stress = inner(f0,Pg_test*f0)
                test_value = project(test_passive_stress-set_point, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
                test_v = interpolate(test_value,Quad)
                test_v_array = test_v.vector().get_local()[:]
                test_max = np.amax(test_v_array)
                print "test max = " +str(test_max)
                #ps = project(test_passive_stress,FunctionSpace(mesh,"DG",1),form_compiler_parameters={"representation":"uflacs"})
                #print "test passive =" +str(ps.vector())"""

        #eff_final = test_E
        #rint "eff final forms " +str(eff_final.vector().get_local()[:])

        #passive stress in fiber direction
        """Sbulk_f = inner(f0,sbulk*f0)
        #passive stress in transverse
        Sbulk_t = inner(s0,sbulk*s0)
        #passive shear stress
        Sbulk_shear = inner(f0,sbulk*s0)

        Sbulk_local = as_tensor([[Sbulk_f, Sbulk_shear, Sbulk_shear],[Sbulk_shear, Sbulk_t, Sbulk_shear],[Sbulk_shear, Sbulk_shear, Sbulk_t]])"""
        #sbulk_flocal = as_tensor([[Sbulk_f, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        #Sbulk_trans_local = as_tensor([[0.0, 0.0, 0.0], [0.0, Sbulk_t, 0.0], [0.0, 0.0, Sbulk_t]])
        #Sbulk_shear_local = as_tensor([[0.0, Sbulk_shear, Sbulk_shear], [Sbulk_shear, 0.0, Sbulk_shear], [Sbulk_shear, Sbulk_shear, 0.0]])
        #S_local = as_tensor([[Sff, Sfs, Sfn], [Sfs, Sss, Sns], [Sfn, Sns, Snn]])

        #S_local = as_tensor([[Sff, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        S_local = as_tensor([[Sff_weighted, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

        TransMatrix = as_tensor(f0[i]*e1[j], (i,j)) + as_tensor(s0[i]*e2[j], (i,j)) + as_tensor(n0[i]*e3[j], (i,j))

        S_global = as_tensor(TransMatrix[i,k]*TransMatrix[j,l]*S_local[k,l],(i,j))

        # guccione global
        """S_bulk_global = as_tensor(TransMatrix[i,k]*TransMatrix[j,l]*Sbulk_local[k,l],(i,j))
        Pg = F*S_global"""
        Pg = F*sbulk
        S = S_global
        #S = S_local

        #P = F*S - p*inv(F.T)
        P = F*S
        Pff =  inner(f0,P*f0)

        #Pbulk = F*S_bulk_global
        #Pbulkf = inner(f0,Pbulk*f0)

        #T = F*S*F.T - p*I
        T = F*S*F.T

        #return  P,S,T, alpha
        #return   Pbulkf, Pff, alpha
        return Pg, Pff, alpha

    def return_radial_vec_ratio(self):

        mesh = self.parameters["mesh"]
        s0 = self.parameters["sheet"]
        print s0[0]

        X = SpatialCoordinate(mesh)
        ratio = s0_evaluated.y()/s0_evaluated.x()

        return ratio

    def Umat(self):

        Fmat = self.Fmat()
        F0 = Fmat
        for j in range(15):
            F0 = 0.5* (F0 + inv(F0).T)
        R = F0
        return inv(R)*Fmat

    def kroon_law(self,FunctionSpace,step_size,kappa):

        mesh = self.parameters["mesh"]
        U = self.Umat()
        f0 = self.parameters["fiber"]
        f = U*f0/sqrt(inner(U*f0,U*f0))
        f_adjusted = 1./kappa * (f - f0) * step_size
        f_adjusted = project(f_adjusted,VectorFunctionSpace(mesh,"DG",1),form_compiler_parameters={"representation":"uflacs"})
        f_adjusted = interpolate(f_adjusted,FunctionSpace)

        return f_adjusted

    def rand_walk(self,width):

        f0 = self.parameters["fiber"]
        for i in np.arange(np.shape(f0.vector().array())[0]/3):
            i = int(i)
            f0.vector()[i*3] = np.random.normal(f0.vector().array()[i*3],width)
            f0.vector()[i*3+1] = np.random.normal(f0.vector().array()[i*3+1],width)
            f0.vector()[i*3+2] = np.random.normal(f0.vector().array()[i*3+2],width)

        return f0
