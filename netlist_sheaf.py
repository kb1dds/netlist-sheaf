import numpy as np
import pysheaf as ps

class NetlistSheaf(ps.Sheaf):
    def __init__(self, parts, nets):
        '''
Construct a netlist sheaf using a dictionary of `parts` and `nets`.  

Each part is named, and has fields:
* `data_dimension`, 
* `bounds` (optional), 
* `ports`.

The `ports` field is a dictionary naming each port, and associating a string that defines a function for the corresponding restriction map.

Example:
```
{ "AND2" : { "data_dimension" : 2,
	     "bounds" : "[(0,1)]*2",
	     "ports" : { "IN1" : "lambda x: x[0]",
		       "IN2" : "lambda x: x[1]",
		       "OUT" : "lambda x: x[0]*x[1]" }},
...}
```

Each net is named, and has fields:
* `data_dimension`,
* `bounds` (optional),
* `connections`.
Subsequently, `connections` is a list of dictionaries, each one has a `part` and `port` field, which index into the `parts` dictionary above.

```
Example
{
    "A": {
	"data_dimension" : 1,
	"bounds" : "[(0,1)]",
	"connections" : [
	       { "part" : "AND2",
		 "port" : "IN1"}
	   ]
    },
...}
```
        '''
        
        ps.Sheaf.__init__(self)

        self.mNumpyNormType=2

        # First, build the parts (lower level of poset)
        for k,v in parts.items():
            self.AddCell(k,
                         ps.Cell('part',
                                 dataDimension = int(v['data_dimension'])))
            self.GetCell(k).SetDataAssignment(ps.Assignment('part',np.zeros((v['data_dimension'],))))

        # Next, build the nets (upper level of poset)
        for k,v in nets.items():
            self.AddCell(k,
                         ps.Cell('net',
                                 dataDimension = int(v['data_dimension'])))
            self.GetCell(k).SetDataAssignment(ps.Assignment('net',np.zeros((v['data_dimension'],))))

            # Connect each part in this net
            for vp in v['connections']:
                part=vp['part']
                port=vp['port']
                self.AddCoface(part,k,
                               ps.Coface('part','net',
                                         eval(parts[part]['ports'][port])))
        return
