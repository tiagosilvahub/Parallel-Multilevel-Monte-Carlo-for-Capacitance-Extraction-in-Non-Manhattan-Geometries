# Parallel-Multilevel-Monte-Carlo-for-Capacitance-Extraction-in-Non-Manhattan-Geometries
A technique for improving the performance of a Monte Carlo algorithm for capacitance extraction of integrated circuits that allowed to reach a x3 better accuracy after an hour of running time. Near perfectly parallel with 99% speedup per CPU core.

Pseudo-code for Walk-On-Spheres algorithm:

```
Input: numberofSimulations, initialP oint, ε, circuit
for numberOfSimulations do
  point = initialPoint
  electrode = NOELECTRODE
  while !collision do
    d, electrode, collision = distance(point, circuit, ε)
    point = particleStep(point, d, collision, electrode)
  end
  v = voltage(electrode, circuit)
  sum = sum + v
end
result = sum/numberOfSimulations
Output: result
```

Pseudo-code for particle step:

```
Input: point, d, collision
rndAngle = rand(0, 2π)
if collision then
  return
end
if electrode == OUT ERLAY ER then
  p.x− = 2 × d × normal.x
  p.y− = 2 × d × normal.y
else
  p.x+ = d × cos(rndAngle)
  p.y+ = d × sin(rndAngle)
end
Output: point
```

Modified Walk-On-Spheres algorithm that has improved performance:
```
Input: numberofSimulations, ε, ε2, level, point, circuit
v1 = v2 = sum = 0
t1 = time()
for numberofSimulations do
  p = initialPoint
  electrode = NOELECTRODE
  first = NOCOLLISION
  second = NOCOLLISION
  while second! = COLLIDED do
    d, first, second, electrode = distance(point, circuit, ε, ε2)
    point = particleStep(point, d, f irst, second, electrode)
    if first == JUSTCOLLIDED then
      if level == 0 then
        second = COLLIDED
        v1 = 0;
      else
        v1 = voltage(electrode, circuit)
        first = COLLIDED
      end
    end
  end
  v2 = voltage(electrode, circuit)
  v = v2 − v1
  sum = sum + v
  sumSq = sumSq + pow(v, 2)
end
t2 = time()
avgCost = (t2 − t1) / numberOfSimulations
avgVariance = sumSq / numberOfSimulations
avgVoltage = sum / numberOfSimulations
Output: avgVoltage, avgVariance, avgCost
```
