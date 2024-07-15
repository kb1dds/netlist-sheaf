import numpy as np
import pysheaf as ps

class NetlistSheaf(ps.Sheaf):
    def __init__(self, parts, nets):
        ps.Sheaf.__init__(self)

        self.mNumpyNormType=2

        for k,v in parts.items():
            self.AddCell(k,
                         ps.Cell('part',
                                 dataDimension = int(v['data_dimension'])))
            self.GetCell(k).SetDataAssignment(ps.Assignment('part',np.zeros((v['data_dimension'],))))

        for k,v in nets.items():
            self.AddCell(k,
                         ps.Cell('net',
                                 dataDimension = int(v['data_dimension'])))
            self.GetCell(k).SetDataAssignment(ps.Assignment('net',np.zeros((v['data_dimension'],))))
            for vp in v['connections']:
                part=vp['part']
                port=vp['port']
                self.AddCoface(part,k,
                               ps.Coface('part','net',
                                         eval(parts[part]['ports'][port])))
        return
