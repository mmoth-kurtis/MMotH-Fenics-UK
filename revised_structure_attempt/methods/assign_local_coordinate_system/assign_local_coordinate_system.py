from dolfin import *
from numpy import random as r

def assign_local_coordinate_system(fiber_fcn_space,lv_options,coord_params):

    # Functions that are useful in unit cube and cylinder for calculating
    # fibrous area, or the long (x) axis to assign local coordinate systems
    test_marker_fcn = Function(marker_space)
    long_axis       = Function(fiberFS)

    if (sim_geometry == "ventricle") or (sim_geometry == "ellipsoid"):
        # we can assign some of these functions from the mesh file
        casename = lv_options["casename"]
        f = lv_options["f"]
        f.read(hsl0, casename + "/" + "hsl0")

        # assign local coordinate system at each gauss point
        f.read(f0, casename + "/" + "eF")
        f.read(s0, casename + "/" + "eS")
        f.read(n0, casename + "/" + "eN")
        f.close()
    else:
        hsl0.vector()[:] = hs_params["initial_hs_length"][0]
        m_x = 1.0
        m_y = 0.0
        m_z = 0.0
        for jj in np.arange(no_of_int_points):

                #if (temp_tester_array[jj] < dig_array[jj]) or (temp_tester_array[jj] > -dig_array[jj] + 10.0):
                if (temp_tester_array[jj] > 9.0) or (temp_tester_array[jj] < 1.0):
                    # inside left end
                    f0.vector()[jj*3] = 1.0
                    f0.vector()[jj*3+1] = 0.0
                    f0.vector()[jj*3+2] = 0.0
                    hs_params_list[jj]["myofilament_parameters"]["k_3"][0] = 0.0
                    #passive_params_list[jj]["c"][0] = 2000
                    c_param.vector()[jj] = 400
                    c2_param.vector()[jj] = 250
                    c3_param.vector()[jj] = 10


                else:

                    # In middle region, assign fiber vector
                    # find radius of point
                    temp_radius = point_rad_array[jj]
                    if np.abs(temp_radius - top_radius) < 0.01:
                        temp_width = 0.0
                    else:
                        temp_width = width*(1.0-(temp_radius*temp_radius/(top_radius*top_radius)))
                    f0.vector()[jj*3] = r.normal(m_x,temp_width,1)[0]
                    f0.vector()[jj*3+1] = r.normal(m_y,temp_width,1)[0]
                    f0.vector()[jj*3+2] = r.normal(m_z,temp_width,1)[0]
                    c_param.vector()[jj] = 1000
                    c2_param.vector()[jj] = 250
                    c3_param.vector()[jj] = 15


    f0 = f0/sqrt(inner(f0,f0))
