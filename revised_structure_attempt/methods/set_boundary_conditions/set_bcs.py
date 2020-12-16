from dolfin import *

## set boundary conditions

def set_bcs(sim_geometry,protocol,mesh,W,facetboundaries,u_D):

    output = {}

    if (sim_geometry == "ventricle") or (sim_geometry == "ellipsoid"):

        # if ventricle or ellipsoid simulation, constrain base in longitudinal direction,
        # and other constraints are in weak form
        topid = 4
        LVendoid = 2
        epiid = 1
        bctop = DirichletBC(W.sub(0).sub(2), Expression(("0.0"), degree = 2), facetboundaries, topid)
        bcs = [bctop]

    elif (sim_geometry == "cylinder") or sim_geometry == "box_mesh" or sim_geometry == "gmesh_cylinder":
        sim_type = protocol["simulation_type"]

        if sim_geometry == "cylinder":
            center = 0.0
            radius = 1.0
        else:
            center = 0.5
            radius = 0.5



        # defining parts of the model where the boundary condition should be applied later
        class Left(SubDomain):
            def inside(self, x, on_boundary):
                tol = 1E-14
                return on_boundary and abs(x[0]) < tol
        #  where x[0] = 10
        class Right(SubDomain):
            def inside(self, x, on_boundary):
                tol = 1E-14
                return on_boundary and abs(x[0]-10.0) < tol
        class Fix_y(SubDomain):
            def inside(self, x, on_boundary):
                tol = 1E-1
                return near(x[0],0.0,tol) and near(x[1],center,tol)
        class Fix_y_right(SubDomain):
            def inside(self, x, on_boundary):
                tol = 1E-14
                return near(x[0],10.0,tol) and near(x[1],center,tol)
        class Fix_z_right(SubDomain):
            def inside(self, x, on_boundary):
                tol = 1E-14
                return near(x[0],10.0,tol) and near(x[2],center,tol)
        class Fix_z(SubDomain):
            def inside(self, x, on_boundary):
                tol = 1E-14
                return (near(x[0],0.0,tol) and near(x[2],center,tol))

        # Appropriately mark all facetboundaries
        facetboundaries.set_all(0)
        left = Left()
        right = Right()
        fix_y = Fix_y()
        fix_y_right = Fix_y_right()
        fix_z = Fix_z()
        fix_z_right = Fix_z_right()

        left.mark(facetboundaries, 1)
        right.mark(facetboundaries, 2)
        fix_y.mark(facetboundaries, 3)
        fix_z.mark(facetboundaries,5)

        # fix left face in x, right face is displacement (until traction bc may be triggered)
        bcleft= DirichletBC(W.sub(0).sub(0), Constant((0.0)), facetboundaries, 1)
        bcright= DirichletBC(W.sub(0).sub(0), u_D, facetboundaries, 2)
        bcfix_y = DirichletBC(W.sub(0).sub(1), Constant((0.0)), fix_y, method="pointwise")
        bcfix_z = DirichletBC(W.sub(0).sub(2), Constant((0.0)), fix_z, method="pointwise")
        bcfix_y_right = DirichletBC(W.sub(0).sub(1), Constant((0.0)),fix_y_right, method="pointwise")
        bcfix_z_right = DirichletBC(W.sub(0).sub(2), Constant((0.0)),fix_z_right, method="pointwise")

        bcs = [bcleft,bcfix_y,bcfix_z,bcfix_y_right,bcfix_z_right,bcright] # order matters!

        if sim_type == "work_loop":
            marker_space = FunctionSpace(mesh,'CG',1)
            bc_right_test = DirichletBC(marker_space,Constant(1),facetboundaries,2)
            test_marker_fcn = Function(marker_space) # this is what we need to grab the displacement after potential shortening
            #bc_right_test.apply(test_marker_fcn.vector())
            output["test_marker_fcn"] = test_marker_fcn

    elif sim_geometry == "unit_cube":
        sim_type = protocol["simulation_type"]

        class Left(SubDomain):
            def inside(self, x, on_boundary):
                tol = 1E-14
                return on_boundary and abs(x[0]) < tol
        #  where x[0] = 10
        class Right(SubDomain):
            def inside(self, x, on_boundary):
                tol = 1E-14
                return on_boundary and abs(x[0]-1.0) < tol
        #  where x[2] = 0
        class Lower(SubDomain):
            def inside(self, x, on_boundary):
                tol = 1E-14
                return on_boundary and abs(x[2]) < tol
        #  where x[1] = 0
        class Front(SubDomain):
            def inside(self, x, on_boundary):
                tol = 1E-14
                return on_boundary and abs(x[1]) < tol
        #  where x[0], x[1] and x[2] = 0
        class Fix(SubDomain):
            def inside(self, x, on_boundary):
                tol = 1E-14
                return on_boundary and abs(x[0]) < tol and abs(x[1]) < tol and abs(x[2]) < tol

        facetboundaries.set_all(0)
        left = Left()
        right = Right()
        fix = Fix()
        lower = Lower()
        front = Front()
        #
        left.mark(facetboundaries, 1)
        right.mark(facetboundaries, 2)
        fix.mark(facetboundaries, 3)
        lower.mark(facetboundaries, 4)
        front.mark(facetboundaries, 5)

        # Similar to cylinder but without fixing displacement along y and z axes to prevent rotation
        bcleft= DirichletBC(W.sub(0).sub(0), Constant((0.0)), facetboundaries, 1)         # u1 = 0 on left face
        bcright= DirichletBC(W.sub(0).sub(0), u_D, facetboundaries, 2)
        bcfix = DirichletBC(W.sub(0), Constant((0.0, 0.0, 0.0)), fix, method="pointwise") # at one vertex u = v = w = 0
        bclower= DirichletBC(W.sub(0).sub(2), Constant((0.0)), facetboundaries, 4)        # u3 = 0 on lower face
        bcfront= DirichletBC(W.sub(0).sub(1), Constant((0.0)), facetboundaries, 5)        # u2 = 0 on front face
        bcs = [bcleft, bclower, bcfront,bcfix, bcright] #order matters!

        if sim_type == "work_loop":
            marker_space = FunctionSpace(mesh,'CG',1)
            bc_right_test = DirichletBC(marker_space,Constant(1),facetboundaries,2)
            test_marker_fcn = Function(marker_space) # this is what we need to grab the displacement after potential shortening
            #bc_right_test.apply(test_marker_fcn.vector())
            output["test_marker_fcn"] = test_marker_fcn

    output["bcs"] = bcs
    return output
