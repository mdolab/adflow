
class DummyOptProb(object):

    def __init__(self):
        self.variables = []
        self.constraints = []
        self.objectives = []

    def addVarGroup(self, *args, **kwargs):
        self.variables.append([args, kwargs])

    def addConGroup(self, *args, **kwargs):
        self.constraints.append([args, kwargs])

    def addVar(self, name, *args, **kwargs):
        self.addVarGroup(name, 1, *args, scalar=True, **kwargs)

    def addCon(self, name, *args, **kwargs):
        self.addConGroup(name, 1, *args, **kwargs)

    def addObj(self, name, *args, **kwargs):
        self.objectives.append([name, args, kwargs])


def get_dvs_and_cons(ap=None, geo=None, con=None):
    vars = []
    cons = []

    if ap is not None:
        ap_vars = DummyOptProb()
        ap.addVariablesPyOpt(ap_vars)
        vars.extend(ap_vars.variables)
        # for long_name, dv in ap.DVs.items(): 
        #     dv_data = {'scalar': True, 
        #                'value': dv.value, 
        #                'lower': dv.lower, 
        #                'upper': dv.upper, 
        #                'scale': dv.scale, 
        #                'offset': dv.offset, 
        #                'units': dv.units} 
        #     vars.append(((dv.key, 1, 'c'), dv_data))

    if geo is not None:
        geo_vars = DummyOptProb()
        geo.addVariablesPyOpt(geo_vars)
        vars.extend(geo_vars.variables)

    if con is not None:
        dv_cons = DummyOptProb()
        con.addConstraintsPyOpt(dv_cons)
        cons.extend(dv_cons.constraints)

    return vars, cons
