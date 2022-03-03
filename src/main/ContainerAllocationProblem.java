package main;
import com.opencsv.CSVReader;
import ec.EvolutionState;
import ec.Individual;
import ec.Population;
import ec.coevolve.GroupedProblemForm;
import ec.gp.GPIndividual;
import ec.gp.GPProblem;
import ec.gp.koza.KozaFitness;
import ec.simple.SimpleEvolutionState;
import ec.util.Parameter;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;


public class ContainerAllocationProblem extends GPProblem implements GroupedProblemForm{
    public static final String P_DATA = "data";
    boolean shouldSetContext;
    MainAllocationProcessContainer mainAllocationProcessContainer;
    MainAllocationProcessVM mainAllocationProcessVM;
    MainAllocationProcessContainer_Ind mainAllocationProcessContainer_Ind;
    MainAllocationProcessContainer_LAH mainAllocationProcessContainer_lah;
    ArrayList<Double> benchmarkResult ;
    SubJustFit_FF benchmark;

    // Initialization data
    ArrayList <ArrayList> initPm;
    ArrayList <ArrayList> initPmType;
    ArrayList <ArrayList> initVm;
    ArrayList <ArrayList> initOs;
    ArrayList <ArrayList> initContainer;

    // A list of containers
    private ArrayList<ArrayList<Double[]>> inputX = new ArrayList<>();

    // A list of candidate VMs, each vm has an array which includes its CPU and Mem capacity
    private ArrayList<Double[]> vmTypeList = new ArrayList<>();
    private ArrayList<Double[]> pmTypeList = new ArrayList<>();

    // An array of OS probability
    private ArrayList<Double> OSPro = new ArrayList<>();

    @Override
    public void setup(final EvolutionState state, final Parameter base){

        // very important
        super.setup(state, base);
        if (!(input instanceof DoubleData)){
            state.output.fatal("GPData class must subclasses from " + DoubleData.class,
                    base.push(P_DATA), null);
        }

        Parameter vmCPUOverheadRateP = new Parameter("VMCPUOverheadRate");
        Parameter vmMemOverheadP = new Parameter("VMMemOverhead");
//        Parameter kP = new Parameter("k");
        Parameter readFileStartFromP = new Parameter("readFileStartFrom");
        Parameter readFileEndP = new Parameter("readFileEnd");
        Parameter testCasePathP = new Parameter("testCasePath");
        Parameter osPathP = new Parameter("osPath");
        Parameter vmConfigPathP = new Parameter("vmConfigPath");
        Parameter pmConfigPathP = new Parameter("pmConfigPath");

        Parameter osProP = new Parameter("osProPath");
        Parameter envPath = new Parameter("initEnvPath");

        double vmCpuOverheadRate = state.parameters.getDouble(vmCPUOverheadRateP, null);
        double vmMemOverhead = state.parameters.getDouble(vmMemOverheadP, null);

//        double k = state.parameters.getDouble(kP, null);
        int start = state.parameters.getInt(readFileStartFromP, null);
        int end = state.parameters.getInt(readFileEndP, null);

        String testCasePath = state.parameters.getString(testCasePathP, null);
        String osPath = state.parameters.getString(osPathP, null);
        String vmConfigPath = state.parameters.getString(vmConfigPathP, null);
        String pmConfigPath = state.parameters.getString(pmConfigPathP, null);
        String osProPath = state.parameters.getString(osProP, null);
        String initEnvPath =  state.parameters.getString(envPath, null);


        readEnvData(initEnvPath, start, end);
        readFromFiles(testCasePath, osPath, start, end - 1);
        readVMConfig(vmConfigPath);
        readPMConfig(pmConfigPath);
        readOSPro(osProPath);


//
//        mainAllocationProcessContainer =
//                new MainAllocationProcessContainer(
//                  this,
//                                        state,
//                                        vmCpuOverheadRate,
//                                        vmMemOverhead
//                                        );

        mainAllocationProcessContainer_Ind =
                new MainAllocationProcessContainer_Ind(
                        this,
                        state,
                        vmCpuOverheadRate,
                        vmMemOverhead
                        );

//        mainAllocationProcessContainer_lah =
//                new MainAllocationProcessContainer_LAH(
//                        this,
//                        state,
//                        vmCpuOverheadRate,
//                        vmMemOverhead
//                );

        mainAllocationProcessVM=
                new MainAllocationProcessVM(
                        this,
                        state,
                        vmCpuOverheadRate,
                        vmMemOverhead
                        );

        benchmark = new SubJustFit_FF(
                vmCpuOverheadRate,
                vmMemOverhead,
                inputX,
                initVm,
                initContainer,
                initOs,
                initPm,
                initPmType,
                pmTypeList,
                vmTypeList
                );

        benchmarkResult = new ArrayList<>();
        for(int i = 0; i < inputX.size(); i++){
            benchmarkResult.add(benchmark.allocate(i));
        }

}


    private void readEnvData(String initEnvPath, int start, int end){
        ReadConfigures readEnvConfig = new ReadConfigures();
        initPm = readEnvConfig.testCases(initEnvPath, "pm", start, end);
        initOs = readEnvConfig.testCases(initEnvPath, "os", start, end);
//        ArrayList<ArrayList> temp= readEnvConfig.testCases(initEnvPath, "container", start, end);
        initContainer = readEnvConfig.testCases(initEnvPath, "container", start, end);
        initVm = readEnvConfig.testCases(initEnvPath, "vm", start, end);
        initPmType = readEnvConfig.testCases(initEnvPath, "pmType", start, end);

    }

    // Read two column from the testCase file
    private ArrayList<Double[]> readFromFile(String path, String osPath){
        ArrayList<Double[]> data = new ArrayList<>();
        try {
            Reader reader = Files.newBufferedReader(Paths.get(path));
            Reader readerOS = Files.newBufferedReader(Paths.get(osPath));

            CSVReader csvReader = new CSVReader(reader);
            CSVReader csvReaderOS = new CSVReader(readerOS);
            String[] nextRecord;
            String[] nextRecordOS;
            while((nextRecord = csvReader.readNext()) != null && (nextRecordOS = csvReaderOS.readNext()) != null){
                // [0] is for cpu, [1] is for mem [2] for os [3] for hours
                Double[] container = new Double[4];
                container[0] = Double.parseDouble(nextRecord[0]);
                container[1] = Double.parseDouble(nextRecord[1]);
                container[2] = Double.parseDouble(nextRecordOS[0]);
                container[3] = Double.parseDouble(nextRecord[2]);
                data.add(container);
            }
            reader.close();
            readerOS.close();
            csvReader.close();
            csvReaderOS.close();
        } catch (IOException e1){
            e1.printStackTrace();
        }
        return data;
    }

    private void readOSPro(String osProPath){
        try {
            Reader reader = Files.newBufferedReader(Paths.get(osProPath));
            CSVReader csvReader = new CSVReader(reader);
            String[] nextRecord;
            while((nextRecord = csvReader.readNext()) != null){
                Double pro;
                pro = Double.parseDouble(nextRecord[0]);
                OSPro.add(pro);
            }
        } catch (IOException e1){
            e1.printStackTrace();
        }
    }

    private void readVMConfig(String vmConfigPath){
        try {
            Reader reader = Files.newBufferedReader(Paths.get(vmConfigPath));
            CSVReader csvReader = new CSVReader(reader);
            String[] nextRecord;
            while((nextRecord = csvReader.readNext()) != null){
                Double[] vm = new Double[3];
                vm[0] = Double.parseDouble(nextRecord[0]); // cpu Mhz
                vm[1] = Double.parseDouble(nextRecord[1]); //mems
                vm[2] = Double.parseDouble(nextRecord[2]); // num of core
                vmTypeList.add(vm);
            }
        } catch (IOException e1){
            e1.printStackTrace();
        }
    }

    private void readPMConfig(String pmConfigPath){
        try {
            Reader reader = Files.newBufferedReader(Paths.get(pmConfigPath));
            CSVReader csvReader = new CSVReader(reader);
            String[] nextRecord;
            while((nextRecord = csvReader.readNext()) != null){
                Double[] pm = new Double[5];
                pm[0] = Double.parseDouble(nextRecord[0]);// cpu
                pm[1] = Double.parseDouble(nextRecord[1]);// mem
                pm[2] = Double.parseDouble(nextRecord[2]);// idle
                pm[3] = Double.parseDouble(nextRecord[3]);// max
                pm[4] = Double.parseDouble(nextRecord[4]);// core num
//                pm[5] = Double.parseDouble(nextRecord[5]); // the type of the pm
                // type will only be used for pm creation
                // due to the creation may not have the full list of the pm types
                // need to this feature to choose the desirable pms
                pmTypeList.add(pm);
            }
        } catch (IOException e1){
            e1.printStackTrace();
        }
    }

    // we read containers from file
    private void readFromFiles(String testCasePath, String osPath, int start, int end){
//        end = 1;
//
        for(int i = start; i <= end; ++i){
//            int toRead = i;
//            if(toRead > 50){
//                toRead = i -50;
//            }
            String path = testCasePath + i + ".csv";
            String pathOS = osPath + i + ".csv";
            inputX.add(readFromFile(path, pathOS));
        }
    }



    @Override
    public void preprocessPopulation(EvolutionState state,
                                     Population pop,
                                     boolean[] prepareForFitnessAssessment,
                                     boolean countVictoriesOnly) {
        for(int i = 0; i < pop.subpops.size(); i++){
            if(prepareForFitnessAssessment[i]){
                for(int j = 0; j < pop.subpops.get(i).individuals.size(); j++){
                    KozaFitness fit = (KozaFitness) (pop.subpops.get(i).individuals.get(j).fitness);
                    fit.trials = new ArrayList();
                }
            }
        }

    }

    @Override
    public int postprocessPopulation(EvolutionState state,
                                     Population pop,
                                     boolean[] assessFitness,
                                     boolean countVictoriesOnly) {
        for(int i = 0; i < pop.subpops.size(); i++){
            if(assessFitness[i]){
                for(int j = 1; j < pop.subpops.get(i).individuals.size(); j++){
                    KozaFitness fit = (KozaFitness)(pop.subpops.get(i).individuals.get(j).fitness);
                    Double fitnessValue = ((Double) fit.trials.get(0)).doubleValue();
                    fit.setStandardizedFitness(state, fitnessValue);
                    pop.subpops.get(i).individuals.get(j).evaluated = true;
                }
            }
        }
        return 0;
    }




    @Override
    public void evaluate(EvolutionState state,
                         Individual[] ind,
                         boolean[] updateFitness,
                         boolean countVictoriesOnly,
                         int[] subpops,
                         int threadnum){
        if (ind.length == 0){
            state.output.fatal("Number of individuals provided to ContainerAllocationProblem is 0");
        }
        if(ind.length == 1){
            state.output.warnOnce("Coevolution used" +
            " but number of individuals provided to RuleCoevolutionProblem is 1.");
        }

        // Step 1: setup VM selection and creation rule
        GPIndividual containerAllocationRule = (GPIndividual) ind[0];

        // Setup VM allocation rule
        GPIndividual vmAllocationRule = (GPIndividual)ind[1];

        KozaFitness fit1 = (KozaFitness) ind[0].fitness;
        KozaFitness fit2 = (KozaFitness) ind[1].fitness;

        if(updateFitness[0]){
//
//            ArrayList<Double> resultList = mainAllocationProcessContainer.evaluate(
//                    (DoubleData)this.input,
//                    state,
//                    inputX,
//                    initVm,
//                    initContainer,
//                    initOs,
//                    initPm,
//                    initPmType,
//                    pmTypeList,
//                    vmTypeList,
//                    containerAllocationRule,
//                    vmAllocationRule,
//                    threadnum, this.stack);

            ArrayList<Double> resultList = mainAllocationProcessContainer_Ind.evaluate(
                    (DoubleData)this.input,
                    state,
                    inputX,
                    initVm,
                    initContainer,
                    initOs,
                    initPm,
                    initPmType,
                    pmTypeList,
                    vmTypeList,
                    containerAllocationRule,
                    vmAllocationRule,
                    threadnum, this.stack);

//            ArrayList<Double> resultList = mainAllocationProcessContainer_lah.evaluate(
//                    (DoubleData)this.input,
//                    state,
//                    inputX,
//                    initVm,
//                    initContainer,
//                    initOs,
//                    initPm,
//                    initPmType,
//                    pmTypeList,
//                    vmTypeList,
//                    containerAllocationRule,
//                    vmAllocationRule,
//                    threadnum, this.stack);

            double aveFit = 0;
            for(int i = 0; i < resultList.size(); ++i){
                aveFit += resultList.get(i);
            }

            // Now we normalize the fitness value
//            aveFit /= benchmarkResult.get(state.generation);
            aveFit /= resultList.size();

            fit1.trials.add(aveFit);
            fit1.setStandardizedFitness(state, aveFit);
//            fit1.setFitness(state, fit1.fitness(), false);
        }

        if(updateFitness[1]){
            ArrayList<Double> resultList = mainAllocationProcessVM.evaluate(
                    (DoubleData)this.input,
                    state,
                    inputX,
                    initVm,
                    initContainer,
                    initOs,
                    initPm,
                    initPmType,
                    pmTypeList,
                    vmTypeList,
                    containerAllocationRule,
                    vmAllocationRule,
                    threadnum, this.stack);

            double aveFit = 0;
            for(int i = 0; i < resultList.size(); ++i){
                aveFit += resultList.get(i);
            }

            // normalize the fitness value
//            aveFit /= benchmarkResult.get(state.generation);
            aveFit /= resultList.size();

            fit2.trials.add(aveFit);
            fit2.setStandardizedFitness(state, aveFit);
//            fit1.setFitness(state, fit1.fitness(), false);
//            System.out.println("Currently we having "+
//                    state.population.subpops.get(0).individuals.size() + " trees ");
//            for(Individual i : state.population.subpops.get(0).individuals) {
//
//            MYGPIndividual myind = (MYGPIndividual) i;
//            ((MYGPIndividual) i).printTrees(state, 1);
//
//        }
        }
    }




    @Override
    public void evaluate(final EvolutionState state,
                         final Individual ind,
                         final int subpopulation,
                         final int threadnum) {
    }

}
