path="AirfoilData.xlsx";
path1="SVSP_velocities.xlsx";
path2="n0012-il.xlsx";
[X,U]=Import_Script(path);
[Utilde,Xtilde]=Import_Script(path1);
CodeCFx(U,X)
figure
CodeCFx(Utilde,abs(Xtilde))
