from dolfin import *
import numpy as np

def update_bcs(bcs,sim_geometry,Ftotal,geo_options,sim_protocol,expr,time,traction_switch_flag,x_dofs):

    output_dict = {}
    print "updating bcs"
    print sim_protocol["simulation_type"]
    # only really need to update if not ventricle simulation
    if (sim_geometry != "ventricle") and (sim_geometry != "ellipsoid"):

        b = assemble(Ftotal,form_compiler_parameters={"representation":"uflacs"})

        for boundary_condition_i in np.arange(np.shape(bcs)[0])-1:
            bcs[boundary_condition_i].apply(b)


        if not geo_options:
            area = 1.0
        elif sim_geometry == "cylinder":
            area = 3.14*geo_options["end_radius"][0]**2 #assuming enough segments are used to approximate a circle
        elif sim_geometry == "box_mesh":
            area = 1.0

        f_int_total = b.copy()
        rxn_force=0.0
        for kk in x_dofs:
            rxn_force += f_int_total[kk]
        output_dict["rxn_force"] = rxn_force

    # If ramp and hold
        # u_D.u_D = whatever it's supposed to be

    if sim_protocol["simulation_type"][0] == "work_loop":
        print "in work loop updating bc"

        if rxn_force >= sim_protocol["traction_magnitude"][0]:
            print "switching to traction boundary condition"
            if traction_switch_flag < 1:
                expr["Press"].P = rxn_force/area
                if sim_geometry == "cylinder" or sim_geometry == "box_mesh":
                    bcs = [bcleft,bcfix_y,bcfix_z,bcfix_y_right,bcfix_z_right]
                if sim_geometry == "unit_cube":
                    bcs = [bcleft, bclower, bcfront,bcfix]
                traction_switch_flag = 1
                output_dict["traction_switch_flag"] = traction_switch_flag
                output_dict["bcs"] = bcs
                output_dict["expressions"]=expr
            else:
                # already switched, keep everything the same
                expr["Press"].P = expr["Press"].P
                bcs = bcs
                output_dict["traction_switch_flag"] = traction_switch_flag
                output_dict["bcs"] = bcs
                output_dict["expr"]=expr
        else:
            # still in the ramp and hold stage
            expr["u_D"].u_D = ramp_and_hold(time,sim_protocol,geo_options)
            output_dict["bcs"] = bcs
            output_dict["expr"]=expr
            output_dict["traction_switch_flag"] = traction_switch_flag

    elif sim_protocol["simulation_type"][0] == "ramp_and_hold":
        expr["u_D"].u_D = ramp_and_hold(time,sim_protocol,geo_options)
        output_dict["expr"] = expr
        output_dict["traction_switch_flag"] = traction_switch_flag
        print "assigning bcs"
        output_dict["bcs"] = bcs

    return output_dict

def ramp_and_hold(time,sim_protocol,geo_options):

    geo_check = not geo_options
    if geo_check:
        # unit cube
        length_scale = 1.0
    else:
        length_scale = geo_options["end_x"][0]
        print "length scale = " + str(length_scale)

    if time < sim_protocol["ramp_t_start"][0]:
        disp = 0.0
        print "displacement is " + str(disp)

    if time >= sim_protocol["ramp_t_end"][0]:
        disp = length_scale*sim_protocol["ramp_magnitude"][0]
        print "displacement is " + str(disp)
    else:
        slope = sim_protocol["ramp_magnitude"][0]/(sim_protocol["ramp_t_end"][0]-sim_protocol["ramp_t_start"][0])
        disp = length_scale*slope*time
        print "displacement is " + str(disp)

    return disp
