# Laminate engineering library for python
A library for python to do laminate calculations. Uses numpy and matplotlib.

## Main features
- Import laminate layup from .csv file
- Plot stress/strain response in x, y, and xy directions
- Possibility to create "Ply" and "Laminate" classes and extract relevant information like stiffness, compliance, thickness ++
- Optinal to set strength and degrading factors for each ply
- Includes a class called "Loadcase" where you can select Laminate and load and check if the current laminate meets the requirements
- By using for loops it is possible to automatically find the best laminate for a given loadcase

##Example

```
plyList = lam.plyListFromCSV('layup.csv')
for ply in plyList:
  ply.setStrength(990, 760, 35, 160, 55)
  ply.setDegradingFactors(1, 0.1, 0.1, 0.5)
lam.plotReactionForce(l, 0,'.\img\plot1')
```
Will create the following image:

![alt tag](https://raw.githubusercontent.com/martinloland/laminate-engineering/master/img/plot1.png)
