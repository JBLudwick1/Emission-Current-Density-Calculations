/*
 * rough_surface_example.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Sep 8 2021, 14:02 by COMSOL 5.6.0.341. */
public class rough_surface_example {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("C:\\Users\\ludwicjb\\Desktop\\MATLAB_code\\Published_LUT_code\\COMSOL_example");

    model.label("rough_surface_example.mph");

    model.param().set("N", "10", "Spatial Frequency Resolution");
    model.param().set("gap", "3", "Anode-Cathode gap");
    model.param().set("Vext", "2.5[kV]", "Anode Voltage");
    model.param().set("L", "1.5e-6", "Isotropic model scaling variable");
    model.param().set("z0", "0.001", "Surface scaling variable");
    model.param().set("phi", "4.5", "Work function in eV");
    model.param().set("temperature", "300 [K]", "Surface temperature");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 3);

    model.result().table().create("tbl1", "Table");

    model.func().create("rn2", "Random");
    model.func().create("extm1", "MATLAB");
    model.func("rn2").label("Uniform Random");
    model.func("rn2").set("funcname", "u1");
    model.func("rn2").set("nargs", 2);
    model.func("rn2").set("uniformrange", "2*pi");
    model.func("extm1").set("clearsolve", true);
    model.func("extm1").set("funcs", new String[][]{{"MATLAB_interpolate", "F,T,R,phi"}});
    model.func("extm1").set("plotargs", new String[][]{{"1", "8"}, {"300", "300"}, {"5", "5"}, {"4", "4"}});

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").geomRep("comsol");
    model.component("comp1").geom("geom1").create("ps2", "ParametricSurface");
    model.component("comp1").geom("geom1").feature("ps2").set("parmax1", 0.5);
    model.component("comp1").geom("geom1").feature("ps2").set("parmax2", 0.5);
    model.component("comp1").geom("geom1").feature("ps2")
         .set("coord", new String[]{"s1*2", "s2*2", "z0*(sum(sum(if((m!=0)||(n!=0),cos(2*pi*(m*s1+n*s2)+u1(m,n)),0),m,-N,N),n,-N,N)+sum(sum(if((m!=0)||(n!=0),cos(2*pi*(m*-s1+n*s2)+u1(m,n)),0),m,-N,N),n,-N,N)+sum(sum(if((m!=0)||(n!=0),cos(2*pi*(m*s1+n*-s2)+u1(m,n)),0),m,-N,N),n,-N,N)+sum(sum(if((m!=0)||(n!=0),cos(2*pi*(m*-s1+n*-s2)+u1(m,n)),0),m,-N,N),n,-N,N))"});
    model.component("comp1").geom("geom1").feature("ps2").set("rtol", "1.0E-3");
    model.component("comp1").geom("geom1").feature("ps2").set("maxknots", 100);
    model.component("comp1").geom("geom1").create("blk2", "Block");
    model.component("comp1").geom("geom1").feature("blk2").set("pos", new String[]{"0", "0", "-gap"});
    model.component("comp1").geom("geom1").feature("blk2").set("size", new String[]{"1", "1", "gap*2"});
    model.component("comp1").geom("geom1").create("pard1", "PartitionDomains");
    model.component("comp1").geom("geom1").feature("pard1").set("partitionwith", "faces");
    model.component("comp1").geom("geom1").feature("pard1").selection("domain").set("blk2(1)", 1);
    model.component("comp1").geom("geom1").feature("pard1").selection("face").set("ps2(1)", 1);
    model.component("comp1").geom("geom1").create("del1", "Delete");
    model.component("comp1").geom("geom1").feature("del1").selection("input").init(3);
    model.component("comp1").geom("geom1").feature("del1").selection("input").set("pard1(1)", 1);
    model.component("comp1").geom("geom1").create("uni1", "Union");
    model.component("comp1").geom("geom1").feature("uni1").selection("input").set("del1", "ps2");
    model.component("comp1").geom("geom1").create("sca1", "Scale");
    model.component("comp1").geom("geom1").feature("sca1").setIndex("factor", "L", 0);
    model.component("comp1").geom("geom1").feature("sca1").selection("input").set("uni1");
    model.component("comp1").geom("geom1").run();
    model.component("comp1").geom("geom1").run("fin");

    model.component("comp1").variable().create("var1");
    model.component("comp1").variable("var1")
         .set("J", "MATLAB_interpolate(es.normE*1e-9,temperature,RoC_adjust,phi)", "local emission current density");
    model.component("comp1").variable("var1")
         .set("RoC", "2e9/(curv1_spatial+curv2_spatial)", "Radius of curvature in nm");
    model.component("comp1").variable("var1")
         .set("RoC_adjust", "if(RoC<0,5e3,RoC)", "Treat concave curvature as a flat surface. This should produce negligible error as emission current from these regions are negligible.");
    model.component("comp1").variable("var1").selection().geom("geom1", 2);
    model.component("comp1").variable("var1").selection().set(2, 3);

    model.component("comp1").material().create("mat1", "Common");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("eta", "Piecewise");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("Cp", "Piecewise");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("rho", "Analytic");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("k", "Piecewise");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("cs", "Analytic");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("an1", "Analytic");
    model.component("comp1").material("mat1").propertyGroup("def").func().create("an2", "Analytic");
    model.component("comp1").material("mat1").propertyGroup().create("RefractiveIndex", "Refractive index");
    model.component("comp1").material("mat1").propertyGroup().create("NonlinearModel", "Nonlinear model");
    model.component("comp1").material("mat1").propertyGroup().create("idealGas", "Ideal gas");
    model.component("comp1").material("mat1").propertyGroup("idealGas").func().create("Cp", "Piecewise");

    model.component("comp1").physics().create("es", "Electrostatics", "geom1");
    model.component("comp1").physics("es").setGroupBySpaceDimension(true);
    model.component("comp1").physics("es").create("gnd1", "Ground", 2);
    model.component("comp1").physics("es").feature("gnd1").selection().set(3);
    model.component("comp1").physics("es").create("pot1", "ElectricPotential", 2);
    model.component("comp1").physics("es").feature("pot1").selection().set(4);

    model.component("comp1").mesh("mesh1").create("ftri1", "FreeTri");
    model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet");
    model.component("comp1").mesh("mesh1").feature("ftri1").selection().set(3);
    model.component("comp1").mesh("mesh1").feature("ftri1").create("size1", "Size");
    model.component("comp1").mesh("mesh1").feature("ftet1").selection().geom("geom1", 3);
    model.component("comp1").mesh("mesh1").feature("ftet1").selection().set(1);
    model.component("comp1").mesh("mesh1").feature("ftet1").create("size1", "Size");

    model.component("comp1").probe().create("bnd2", "Boundary");
    model.component("comp1").probe("bnd2").selection().set(3);

    model.result().table("tbl1").label("Probe Table 1");

    model.component("comp1").material("mat1").label("Air");
    model.component("comp1").material("mat1").set("family", "air");
    model.component("comp1").material("mat1").propertyGroup("def").func("eta").set("arg", "T");
    model.component("comp1").material("mat1").propertyGroup("def").func("eta")
         .set("pieces", new String[][]{{"200.0", "1600.0", "-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4"}});
    model.component("comp1").material("mat1").propertyGroup("def").func("eta").set("argunit", "K");
    model.component("comp1").material("mat1").propertyGroup("def").func("eta").set("fununit", "Pa*s");
    model.component("comp1").material("mat1").propertyGroup("def").func("Cp").set("arg", "T");
    model.component("comp1").material("mat1").propertyGroup("def").func("Cp")
         .set("pieces", new String[][]{{"200.0", "1600.0", "1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4"}});
    model.component("comp1").material("mat1").propertyGroup("def").func("Cp").set("argunit", "K");
    model.component("comp1").material("mat1").propertyGroup("def").func("Cp").set("fununit", "J/(kg*K)");
    model.component("comp1").material("mat1").propertyGroup("def").func("rho")
         .set("expr", "pA*0.02897/R_const[K*mol/J]/T");
    model.component("comp1").material("mat1").propertyGroup("def").func("rho").set("args", new String[]{"pA", "T"});
    model.component("comp1").material("mat1").propertyGroup("def").func("rho").set("argunit", "Pa,K");
    model.component("comp1").material("mat1").propertyGroup("def").func("rho").set("fununit", "kg/m^3");
    model.component("comp1").material("mat1").propertyGroup("def").func("rho")
         .set("plotargs", new String[][]{{"pA", "0", "1"}, {"T", "0", "1"}});
    model.component("comp1").material("mat1").propertyGroup("def").func("k").set("arg", "T");
    model.component("comp1").material("mat1").propertyGroup("def").func("k")
         .set("pieces", new String[][]{{"200.0", "1600.0", "-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4"}});
    model.component("comp1").material("mat1").propertyGroup("def").func("k").set("argunit", "K");
    model.component("comp1").material("mat1").propertyGroup("def").func("k").set("fununit", "W/(m*K)");
    model.component("comp1").material("mat1").propertyGroup("def").func("cs")
         .set("expr", "sqrt(1.4*R_const[K*mol/J]/0.02897*T)");
    model.component("comp1").material("mat1").propertyGroup("def").func("cs").set("args", new String[]{"T"});
    model.component("comp1").material("mat1").propertyGroup("def").func("cs").set("argunit", "K");
    model.component("comp1").material("mat1").propertyGroup("def").func("cs").set("fununit", "m/s");
    model.component("comp1").material("mat1").propertyGroup("def").func("cs")
         .set("plotargs", new String[][]{{"T", "273.15", "373.15"}});
    model.component("comp1").material("mat1").propertyGroup("def").func("an1").set("funcname", "alpha_p");
    model.component("comp1").material("mat1").propertyGroup("def").func("an1")
         .set("expr", "-1/rho(pA,T)*d(rho(pA,T),T)");
    model.component("comp1").material("mat1").propertyGroup("def").func("an1").set("args", new String[]{"pA", "T"});
    model.component("comp1").material("mat1").propertyGroup("def").func("an1").set("argunit", "Pa,K");
    model.component("comp1").material("mat1").propertyGroup("def").func("an1").set("fununit", "1/K");
    model.component("comp1").material("mat1").propertyGroup("def").func("an1")
         .set("plotargs", new String[][]{{"pA", "101325", "101325"}, {"T", "273.15", "373.15"}});
    model.component("comp1").material("mat1").propertyGroup("def").func("an2").set("funcname", "muB");
    model.component("comp1").material("mat1").propertyGroup("def").func("an2").set("expr", "0.6*eta(T)");
    model.component("comp1").material("mat1").propertyGroup("def").func("an2").set("args", new String[]{"T"});
    model.component("comp1").material("mat1").propertyGroup("def").func("an2").set("argunit", "K");
    model.component("comp1").material("mat1").propertyGroup("def").func("an2").set("fununit", "Pa*s");
    model.component("comp1").material("mat1").propertyGroup("def").func("an2")
         .set("plotargs", new String[][]{{"T", "200", "1600"}});
    model.component("comp1").material("mat1").propertyGroup("def").set("thermalexpansioncoefficient", "");
    model.component("comp1").material("mat1").propertyGroup("def").set("molarmass", "");
    model.component("comp1").material("mat1").propertyGroup("def").set("bulkviscosity", "");
    model.component("comp1").material("mat1").propertyGroup("def")
         .set("thermalexpansioncoefficient", new String[]{"alpha_p(pA,T)", "0", "0", "0", "alpha_p(pA,T)", "0", "0", "0", "alpha_p(pA,T)"});
    model.component("comp1").material("mat1").propertyGroup("def").set("molarmass", "0.02897[kg/mol]");
    model.component("comp1").material("mat1").propertyGroup("def").set("bulkviscosity", "muB(T)");
    model.component("comp1").material("mat1").propertyGroup("def")
         .set("relpermeability", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    model.component("comp1").material("mat1").propertyGroup("def")
         .set("relpermittivity", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    model.component("comp1").material("mat1").propertyGroup("def").set("dynamicviscosity", "eta(T)");
    model.component("comp1").material("mat1").propertyGroup("def").set("ratioofspecificheat", "1.4");
    model.component("comp1").material("mat1").propertyGroup("def")
         .set("electricconductivity", new String[]{"0[S/m]", "0", "0", "0", "0[S/m]", "0", "0", "0", "0[S/m]"});
    model.component("comp1").material("mat1").propertyGroup("def").set("heatcapacity", "Cp(T)");
    model.component("comp1").material("mat1").propertyGroup("def").set("density", "rho(pA,T)");
    model.component("comp1").material("mat1").propertyGroup("def")
         .set("thermalconductivity", new String[]{"k(T)", "0", "0", "0", "k(T)", "0", "0", "0", "k(T)"});
    model.component("comp1").material("mat1").propertyGroup("def").set("soundspeed", "cs(T)");
    model.component("comp1").material("mat1").propertyGroup("def").addInput("temperature");
    model.component("comp1").material("mat1").propertyGroup("def").addInput("pressure");
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("n", "");
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("ki", "");
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("n", "");
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("ki", "");
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("n", "");
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("ki", "");
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex")
         .set("n", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex")
         .set("ki", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    model.component("comp1").material("mat1").propertyGroup("NonlinearModel").set("BA", "(def.gamma+1)/2");
    model.component("comp1").material("mat1").propertyGroup("idealGas").func("Cp").label("Piecewise 2");
    model.component("comp1").material("mat1").propertyGroup("idealGas").func("Cp").set("arg", "T");
    model.component("comp1").material("mat1").propertyGroup("idealGas").func("Cp")
         .set("pieces", new String[][]{{"200.0", "1600.0", "1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4"}});
    model.component("comp1").material("mat1").propertyGroup("idealGas").func("Cp").set("argunit", "K");
    model.component("comp1").material("mat1").propertyGroup("idealGas").func("Cp").set("fununit", "J/(kg*K)");
    model.component("comp1").material("mat1").propertyGroup("idealGas").set("Rs", "R_const/Mn");
    model.component("comp1").material("mat1").propertyGroup("idealGas").set("heatcapacity", "Cp(T)");
    model.component("comp1").material("mat1").propertyGroup("idealGas").set("ratioofspecificheat", "1.4");
    model.component("comp1").material("mat1").propertyGroup("idealGas").set("molarmass", "0.02897");
    model.component("comp1").material("mat1").propertyGroup("idealGas").addInput("temperature");
    model.component("comp1").material("mat1").propertyGroup("idealGas").addInput("pressure");

    model.component("comp1").physics("es").feature("pot1").set("V0", "Vext");

    model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hauto", 1);
    model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hmax", "L/50");
    model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hmaxactive", true);
    model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hmin", 6.38E-11);
    model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hminactive", false);
    model.component("comp1").mesh("mesh1").run();

    model.component("comp1").probe("bnd2").label("Emisison Current (A)");
    model.component("comp1").probe("bnd2").set("type", "integral");
    model.component("comp1").probe("bnd2").set("probename", "Ie");
    model.component("comp1").probe("bnd2").set("expr", "J");
    model.component("comp1").probe("bnd2").set("unit", "");
    model.component("comp1").probe("bnd2").set("descr", "");
    model.component("comp1").probe("bnd2").set("table", "tbl1");
    model.component("comp1").probe("bnd2").set("window", "window1");

    model.study().create("std1");
    model.study("std1").create("stat", "Stationary");

    model.sol().create("sol1");
    model.sol("sol1").study("std1");
    model.sol("sol1").attach("std1");
    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol1").create("s1", "Stationary");
    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").create("i1", "Iterative");
    model.sol("sol1").feature("s1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol1").feature("s1").feature().remove("fcDef");
    model.sol().create("sol2");
    model.sol("sol2").study("std1");
    model.sol("sol2").label("Parametric Solutions 1");

    model.result().dataset().create("surf1", "Surface");
    model.result().dataset().create("int2", "Integral");
    model.result().dataset().create("dset3", "Solution");
    model.result().dataset("surf1").selection().set(3);
    model.result().dataset("dset2").set("probetag", "bnd2");
    model.result().dataset("dset2").set("solution", "sol1");
    model.result().dataset("int2").set("probetag", "bnd2");
    model.result().dataset("int2").set("data", "dset2");
    model.result().dataset("int2").selection().geom("geom1", 2);
    model.result().dataset("int2").selection().set(3);
    model.result().dataset("dset3").set("solution", "sol2");
    model.result().numerical().create("pev2", "EvalPoint");
    model.result().numerical("pev2").set("probetag", "bnd2");
    model.result().create("pg4", "PlotGroup3D");
    model.result().create("pg5", "PlotGroup1D");
    model.result("pg4").create("surf1", "Surface");
    model.result("pg4").feature("surf1").set("expr", "if(log10(J)<-50,-50,log10(J))");
    model.result("pg5").set("probetag", "window1");
    model.result("pg5").create("tblp1", "Table");
    model.result("pg5").feature("tblp1").set("probetag", "bnd2");

    model.component("comp1").probe("bnd2").genResult(null);

    model.result("pg5").tag("pg5");

    model.sol("sol1").attach("std1");
    model.sol("sol1").feature("st1").label("Compile Equations: Stationary");
    model.sol("sol1").feature("v1").label("Dependent Variables 1.1");
    model.sol("sol1").feature("s1").label("Stationary Solver 1.1");
    model.sol("sol1").feature("s1").feature("dDef").label("Direct 1");
    model.sol("sol1").feature("s1").feature("aDef").label("Advanced 1");
    model.sol("sol1").feature("s1").feature("fc1").label("Fully Coupled 1.1");
    model.sol("sol1").feature("s1").feature("i1").label("Iterative 1.1");
    model.sol("sol1").feature("s1").feature("i1").set("linsolver", "cg");
    model.sol("sol1").feature("s1").feature("i1").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").label("Multigrid 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("prefun", "amg");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("dDef").label("Direct 1");
    model.sol("sol1").runAll();

    model.result().dataset("dset2").label("Probe Solution 2");
    model.result().numerical("pev2").set("descr", new String[]{""});
    model.result().numerical("pev2").setResult();
    model.result("pg4").label("Local Emission Current Density Plot");
    model.result("pg4").feature("surf1").label("Local Emission Current Density Plot");
    model.result("pg4").feature("surf1").set("data", "surf1");
    model.result("pg4").feature("surf1").set("descractive", true);
    model.result("pg4").feature("surf1").set("descr", "J in logspace. All 0 currents converted to log10(J)=-50");
    model.result("pg4").feature("surf1").set("resolution", "normal");
    model.result("pg5").set("xlabel", "J, Emisison Current (A)");
    model.result("pg5").set("ylabel", "J, Emisison Current (A)");
    model.result("pg5").set("xlabelactive", false);
    model.result("pg5").set("ylabelactive", false);

    return model;
  }

  public static void main(String[] args) {
    run();
  }

}
