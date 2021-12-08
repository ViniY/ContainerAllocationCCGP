package main;

import com.opencsv.CSVReader;
import ec.EvolutionState;
import ec.Individual;
import ec.gp.ADFStack;
import ec.gp.GPIndividual;
import ec.simple.SimpleEvolutionState;
//import sun.java2d.pipe.SpanShapeRenderer;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLOutput;
import java.util.ArrayList;
import java.util.HashMap;

public class MainAllocationProcessContainer {
    public double maxCPU = 211200;
    public double maxMem = 4096000;
    double maxVMMem = 128000;
    public final double PMENERGY;



    private final double vmCpuOverheadRate;
    private final double vmMemOverhead;

    public double containerCpu;
    public double containerMem;
    public int containerOs;
    public double containerOsPro;


    public double normalizedContainerCpu;
    public double normalizedContainerMem;
    public double normalizedVmCpuCapacity;
    public double normalizedVmMemCapacity;


    public double currentPmCpuRemain;
    public double currentPmMemRemain;

    public double normalizedVmCpuOverhead;
    public double normalizedVmMemOverhead;

    public double normalizedGlobalCpuWaste;
    public double normalizedGlobalMemWaste;

    public double normalizedPmCpuRemain;
    public double normalizedPmMemRemain;

    public double normalizedPmActualCpuUsed;
    public double normalizedPmActualMemUsed;

    public double normalizedVmActualCpuUsed;
    public double normalizedVmActualMemUsed;



    public boolean newVmFlag = false;


    private final double k;


    private ContainerAllocationProblem containerAllocationProblem;

    public ArrayList<Integer> rentedPM = new ArrayList<>();
    // Constructor
    public MainAllocationProcessContainer(
                                    ContainerAllocationProblem containerAllocationProblem,
                                    EvolutionState state,
                                    double PMENERGY,
                                    double vmCpuOverheadRate,
                                    double vmMemOverhead,
                                    double k)
    {

        this.containerAllocationProblem = containerAllocationProblem;
        this.PMENERGY = PMENERGY;
        this.vmCpuOverheadRate = vmCpuOverheadRate;
        this.vmMemOverhead = vmMemOverhead;
        this.k = k;
        this.rentedPM = new ArrayList<>();

        // Initialize state


//        MyEvolutionState myEvolutionState = new MyEvolutionState();

        MyEvolutionState myEvolutionState = (MyEvolutionState) state;

//        myEvolutionState.1.0 = 1.0;
//        myEvolutionState.1.0 = 1.0;
        myEvolutionState.PMENERGY = PMENERGY;

    }


    private void initializeDataCenter(int testCase,
                                      ArrayList<Double[]> pmResourceList,
                                      ArrayList<Double[]> pmActualUsageList,
                                      ArrayList<Double[]> vmResourceList,
                                      ArrayList<Double[]> vmTypeList,
                                      ArrayList<Double[]> pmTypeList,
                                      ArrayList<ArrayList> initPm,
                                      ArrayList<ArrayList> initVm,
                                      ArrayList<ArrayList> initContainer,
                                      ArrayList<ArrayList> initOs,
                                      ArrayList<ArrayList> initPmType,
                                      HashMap<Integer, Integer> VMPMMapping,
                                      HashMap<Integer, Integer> vmIndexTypeMapping){

        ArrayList<Double[]> initPmList = initPm.get(testCase);
        ArrayList<Double[]> initVmList = initVm.get(testCase);
        ArrayList<Double[]> containerList = initContainer.get(testCase);
        ArrayList<Double[]> osList = initOs.get(testCase);
        ArrayList<Double[]> initPmTypeList = initPmType.get(testCase);

        int globalVmCounter = 0;
        for(int i =0; i < initPmTypeList.size(); ++i){
            int typePM =(initPmTypeList.get(i)[0]).intValue();
            Double[] vms = initPmList.get(i);
           double pmCPU = pmTypeList.get(typePM)[0];
           double pmMem = pmTypeList.get(typePM)[1];
           // Add the PM to resource List at the beginning stage
            // The order of elements in the list is CPU, Memory and the type of this PM
           pmResourceList.add(new Double[]{pmCPU,pmMem,(Double)(typePM*1.0)});
           pmActualUsageList.add(new Double[]{pmCPU,pmMem});
           // for this vm
            for (int vmCounter = 0; vmCounter <vms.length; ++vmCounter ){
                // Get the type of this VM
                int vmType = vms[vmCounter].intValue() ;

                // Get the OS type
                Double[] os = osList.get(vmCounter + globalVmCounter);

                // Create this VM
                vmResourceList.add(new Double[]{
                        vmTypeList.get(vmType)[0] - vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
                        vmTypeList.get(vmType)[1] - vmMemOverhead,
                        new Double(os[0])
                });
                // get the containers allocated on this VM
                Double[] containers = initVmList.get(vmCounter + globalVmCounter);

                // Allocate the VM to this PM,
                // Allocation includes two part, first, pmResourceList indicates the left resource of PM (subtract entire VMs' size)
                // pmIndex denotes the last PM. pmIndex should be at least 0.
                int pmIndex = pmResourceList.size() - 1;

                // update the pm left resources
                pmResourceList.set(pmIndex, new Double[]{
                        pmResourceList.get(pmIndex)[0] - vmTypeList.get(vmType)[0],
                        pmResourceList.get(pmIndex)[1] - vmTypeList.get(vmType)[1]
                });

                // The second part of allocation,
                // We update the actual usage of PM's resources
                pmActualUsageList.set(pmIndex, new Double[]{
                        pmActualUsageList.get(pmIndex)[0] - vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
                        pmActualUsageList.get(pmIndex)[1] - vmMemOverhead
                });

                // Map the VM to the PM
                VMPMMapping.put(vmCounter + globalVmCounter, pmIndex);

                // for each container
                for(int conContainer = containers[0].intValue();
                    conContainer < containers[containers.length - 1].intValue();
                    ++conContainer){

                    // Get the container's cpu and memory
                    Double[] cpuMem = containerList.get(conContainer);

                    //Create this container
                    // get the left resources of this VM
                    int vmIndex = vmResourceList.size() - 1;
                    Double[] vmCpuMem = vmResourceList.get(vmIndex);

                    // update the vm
                    vmResourceList.set(vmIndex, new Double[] {
                            vmCpuMem[0] - cpuMem[0],
                            vmCpuMem[1] - cpuMem[1],
                            new Double(os[0])
                    });
                    // Whenever we create a new VM, map its index in the VMResourceList to its type for future purpose
                    for (int j = 0; j < vmResourceList.size(); j++) {
                        vmIndexTypeMapping.put(j, vmResourceList.get(j)[2].intValue());
                    }
//                    vmIndexTypeMapping.put(vmResourceList.size() - 1, vmType);
                    // Add the Actual usage to the PM
                    // Here, we must consider the overhead
                    Double[] pmCpuMem = pmActualUsageList.get(pmIndex);

                    // update the pm
                    pmActualUsageList.set(pmIndex, new Double[]{
                            pmCpuMem[0] - cpuMem[0],
                            pmCpuMem[1] - cpuMem[1]
                    });


//                    vm.addContainer(container);
                } // Finish allocate containers to VMs
            } // End  of each VM
            // we must update the globalVmCounter
            globalVmCounter += vms.length;

        }// End  of each PM

//        System.out.println(globalVmCounter);
    }




    /**
     * Calculate the energy consumption using the following equation:
     * Energy = k * MaxEnergy + (1 - k)  * MaxEnergy * utilization_of_a_PM
     * @param pmActualUsageList
     * @return the energy consumption
     */
    private Double energyCalculation(ArrayList<Double[]> pmActualUsageList, ArrayList<Double[]> pmResourceList){
        Double energy = 0.0;
        for (int i = 0; i < pmResourceList.size(); i++) {
            energy+= (pmActualUsageList.get(i)[0])/maxCPU *PMENERGY *(1-k) + k*PMENERGY * maxCPU;
        }
//        for(Double[] pmActualResource:pmActualUsageList){
//            energy += ((PMCPU - pmActualResource[0]) / 1.0) * PMENERGY * (1 - k) + k * PMENERGY * this.maxCPU;
//        energy += ((PMCPU - pmActualResource[0]) / PMCPU) * PMENERGY * (1 - k) + k * PMENERGY;
//        }
        return energy;
    }
//
    private Double pmUtilRemain(ArrayList<Double[]> pmActualUsageList, ArrayList<Double[]> vmResourceList){
        Double meanPmUtil = 0.0;
        for(Double[] pmActualResource:pmActualUsageList){
            meanPmUtil += pmActualResource[0] / 1.0;
        }

        Double aveMemOverhead = vmResourceList.size() * vmMemOverhead / pmActualUsageList.size() / 1.0;
//        Double meanPmResource = 0.0;
//        for(Double[] pmResource:pmResourceList){
//            meanPmResource += (pmResource[0] / 1.0);
//        }
        meanPmUtil /= pmActualUsageList.size();
//        meanPmResource /= pmResourceList.size();
//        System.out.println("meanPmUtil = " + meanPmUtil);
//        System.out.println("aveMemOverhead = " + aveMemOverhead);
//
//
//        System.out.println("meanPmUtil: " + meanPmUtil);
//        System.out.println("aveMemOverhead: " + aveMemOverhead);
//        System.out.println("vmResourceList.size() = " + vmResourceList.size() + " pmActualUsageList.size() = " + pmActualUsageList.size());
        return meanPmUtil + aveMemOverhead;
    }
//    private Double vmUtil(ArrayList<Double[]> vmResourceList){
//        Double meanVmUtil = 0.0;
//        for(Double[] vmResource:vmResourceList){
//            meanVmUtil += vmResource[0];
//        }
//        meanVmUtil /= vmResourceList.size();
//        return meanVmUtil;
//    }

    public ArrayList<Double> evaluate(
                                DoubleData input,
                                EvolutionState state,
                                ArrayList<ArrayList<Double[]>> inputX,
                                ArrayList<ArrayList> initVm,
                                ArrayList<ArrayList> initContainer,
                                ArrayList<ArrayList> initOs,
                                ArrayList<ArrayList> initPm,
                                ArrayList<ArrayList> initPmType,
                                ArrayList<Double[]> pmTypeList,
                                ArrayList<Double[]> vmTypeList,
                                GPIndividual vmSelectionCreationRule,
                                GPIndividual vmAllocationRule,
                                int threadnum,
                                ADFStack stack){

        MyEvolutionState myEvolutionState = (MyEvolutionState) state;
        // testCaseNum equals the current generation
        int testCase = state.generation;
        // initialize the resource lists
        ArrayList<Double> resultList = new ArrayList<>();
        ArrayList<Double> comparedResultList = new ArrayList<>();

        // Loop through the testCases
//        for (int testCase = 0; testCase <= end - start - 1; ++testCase) {


            double globalCPUWaste = 0;
            double globalMEMWaste = 0;
            ArrayList<Double[]> pmResourceList = new ArrayList<>();
            ArrayList<Double[]> pmActualUsageList = new ArrayList<>();
            ArrayList<Double[]> vmResourceList = new ArrayList<>();
            HashMap<Integer, Integer> VMPMMapping = new HashMap<>();
            HashMap<Integer, Integer> vmIndexTypeMapping = new HashMap<>();

            // Initialize data center
            initializeDataCenter(
                    testCase,
                    pmResourceList,
                    pmActualUsageList,
                    vmResourceList,
                    vmTypeList,
                    pmTypeList,
                    initPm,
                    initVm,
                    initContainer,
                    initOs,
                    initPmType,
                    VMPMMapping,
                    vmIndexTypeMapping);

            // the total energy
            Double Energy = energyCalculation(pmActualUsageList, pmResourceList);
//            Double pmFit = pmUtilRemain(pmActualUsageList, pmResourceList);


            // No data center initialization
//                Double Energy = 0.0;
            // Get all the containers
            ArrayList<Double[]> containers = inputX.get(testCase);


            int containerCounter = 0;
            // Start simulation
            for (Double[] container:containers) {

                containerCpu = container[0];
                containerMem = container[1];
                containerOs = container[2].intValue();

                // update myEvolutionState
                myEvolutionState.containerCpu = containerCpu;
                myEvolutionState.containerMem = containerMem;
                myEvolutionState.containerOs = containerOs;




                Integer chosenVM;
                Integer currentVmNum = vmResourceList.size();
                for (int i = 0; i < vmResourceList.size(); i++) {
                    vmIndexTypeMapping.put(i, vmResourceList.get(i)[2].intValue());
                }

                // select or create a VM
                chosenVM = VMSelectionCreation(
                        input,
                        state,
                        vmSelectionCreationRule,
                        threadnum,
                        stack,
                        vmResourceList,
                        pmResourceList,
                        pmActualUsageList,
                        vmTypeList,
                        vmIndexTypeMapping,
                        containerCpu,
                        containerMem,
                        containerOs,
                        globalCPUWaste,
                        globalMEMWaste);

                // check if the VM exists, if chosenVM < currentVmNum is true, it means
                // the chosenVM exists, we just need to update its resources

                if (chosenVM < currentVmNum) {
                    // update the VM resources, allocating this container into this VM
                    vmResourceList.set(chosenVM, new Double[]{
                            vmResourceList.get(chosenVM)[0] - containerCpu,
                            vmResourceList.get(chosenVM)[1] - containerMem,
                            new Double(containerOs)
                    });

                    // Find the pmIndex in the mapping
                    int pmIndex = VMPMMapping.get(chosenVM);

                    // update the PM actual resources
                    pmActualUsageList.set(pmIndex, new Double[]{
                            pmActualUsageList.get(pmIndex)[0] - containerCpu,
                            pmActualUsageList.get(pmIndex)[1] - containerMem
                    });

                    // Else, we need to create this new VM
                } else {

                    // Retrieve the type of select VM
                    int vmType = chosenVM - currentVmNum;

                    // create this new VM
                    vmResourceList.add(new Double[]{
                            vmTypeList.get(vmType)[0] - containerCpu - vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
                            vmTypeList.get(vmType)[1] - containerMem - vmMemOverhead,
                            new Double(containerOs)
                    });

                    // Whenever we create a new VM, map its index in the VMResourceList to its type for future purpose
                    vmIndexTypeMapping = new HashMap<>();
                    for (int i = 0; i < vmResourceList.size() ; i++) {
                        vmIndexTypeMapping.put(i,vmResourceList.get(i)[2].intValue());
                    }
//                    vmIndexTypeMapping.put(vmResourceList.size() - 1, vmType);

                    // After creating a VM, we will choose a PM to allocate
                    Integer chosenPM = VMAllocation(
                                            input,
                                            myEvolutionState,
                                            vmAllocationRule,
                                            threadnum,
                                            stack,
                                            pmResourceList,
                                            pmActualUsageList,
                                            pmTypeList,
                                            vmTypeList.get(vmType)[0],
                                            vmTypeList.get(vmType)[1],
                                  containerCpu + vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
                                  containerMem + vmMemOverhead
                                            );

                    // If we cannot choose a PM
                    if (chosenPM == null) {

                        // Add the VM to the newly created PM
                        // We don't need to consider the overhead here.
                        Double[] chosedPM = pmCreation(pmTypeList,vmTypeList.get(vmType)[0],vmTypeList.get(vmType)[1] );
                        pmResourceList.add(chosedPM);

                        // Add the Actual usage to the PM
                        // Here, we must consider the overhead
                        pmActualUsageList.add(new Double[]{
                                chosedPM[0] - containerCpu - vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
                                chosedPM[1] - containerMem - vmMemOverhead
                        });

                        // Map the VM to the PM
                        VMPMMapping.put(vmResourceList.size() - 1, pmResourceList.size() - 1);


                        // If there is an existing PM, we allocate it to an existing PM
                    } else {

                        currentPmCpuRemain = pmResourceList.get(chosenPM)[0] - vmTypeList.get(vmType)[0];
                        currentPmMemRemain = pmResourceList.get(chosenPM)[1] - vmTypeList.get(vmType)[1];

                        // update the PM resources
                        // pm resources - vm size
                        pmResourceList.set(chosenPM, new Double[]{
                                currentPmCpuRemain,
                                currentPmMemRemain
                        });

                        // update the actual resources
                        // Actual usage - container required - vm overhead
                        pmActualUsageList.set(chosenPM, new Double[]{
                                pmActualUsageList.get(chosenPM)[0] - containerCpu - vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
                                pmActualUsageList.get(chosenPM)[1] - containerMem - vmMemOverhead
                        });

                        // Map the VM to the PM
                        VMPMMapping.put(vmResourceList.size() - 1, chosenPM);

                    } // End of allocating a VM to an existing PM

                } // End of creating a new VM



                Energy = energyCalculation(pmActualUsageList, pmResourceList);

                resultList.add(Energy);
        } // End of all test cases


        return resultList;
    }
    private Double[] pmCreation(ArrayList<Double[]> pmTypeList, double vmCpu, double vmMem) {
        int chosedType = -1;
        double bestFitValue = 0;
        double bestcurrentUtil_CPU =0;
        double bestcurrentUtil_Mem = 0;
        double requireCPU = vmCpu;
        double requireMem = vmMem;

        for (int i = 0; i < pmTypeList.size(); i++) {
            if(requireCPU < pmTypeList.get(i)[0] && requireMem < pmTypeList.get(i)[1]) {
                double currentUtil_CPU = (pmTypeList.get(i)[0] -requireCPU ) / pmTypeList.get(i)[0];
                double currentUtil_Mem = (pmTypeList.get(i)[1] - requireMem) / pmTypeList.get(i)[1];
                double fitValue =  bestcurrentUtil_CPU * bestcurrentUtil_Mem;
                if (bestcurrentUtil_CPU < currentUtil_CPU) {
                    chosedType = i;
                    bestcurrentUtil_CPU = currentUtil_CPU;
//                    bestcurrentUtil_Mem = currentUtil_Mem;
//                    bestFitValue = fitValue;
                }
//                else if (fitValue > bestFitValue){// do not have both requirements better than others but together is better
//                    bestcurrentUtil_Mem  = currentUtil_Mem;
//                    bestcurrentUtil_CPU = currentUtil_CPU;
//                    bestFitValue = fitValue;
//                    chosedType = i;
//                }
                else{//the current type is worse
                    continue;
                }
            }

        }
        if (chosedType <0 || chosedType>pmTypeList.size()) chosedType =0;
        Double newPM[] = new Double[]{pmTypeList.get(chosedType)[0], pmTypeList.get(chosedType)[1], chosedType*1.0};
        return newPM;
    }



//    /**
//     * VMAllocation allocates a VM to a PM
//     * @param pmResourceList is a list of remaining resources of PM, not the acutal used resources
//     * @param vmCpuCapacity the cpu capacity of a VM
//     * @param vmMemCapacity the memory capacity of a VM
//     * @return an index of a PM
//     */
//    private Integer VMAllocation(ArrayList<Double[]> pmResourceList, double vmCpuCapacity, double vmMemCapacity){
//        Integer chosenPM = null;
////        Double BestScore = null;
//
//        // Loop through the PMs in the existing PM list
//        for(int pmCount = 0; pmCount < pmResourceList.size(); ++pmCount){
//            double pmCpuRemain = pmResourceList.get(pmCount)[0];
//            double pmMemRemain = pmResourceList.get(pmCount)[1];
//
//            // First Fit
//            if (pmCpuRemain >= vmCpuCapacity && pmMemRemain >= vmMemCapacity) {
//                chosenPM = pmCount;
//                break;
//            } // End if
//        }
//
//        return chosenPM;
//    }


    /*
     * Basically, VMAllocation still use the BestFit framework
     */
    private Integer VMAllocation(
            DoubleData input,
            final EvolutionState state,
            final Individual ind,
            final int threadnum,
            final ADFStack stack,
            ArrayList<Double[]> pmResourceList,
            ArrayList<Double[]> pmActualResourceList,
            ArrayList<Double[]> pmTypeList,
            double vmCpuCapcacity,
            double vmMemCapacity,
            double vmUsedCpu,
            double vmUsedMem
            ){

        Integer chosenPM = null;
        Double BestScore = null;
        int pmCount = 0;
        // Loop through the tempResourceList
        for(Double[] pm:pmResourceList){
            // Get the remaining PM resources
            double pmCpuRemain = pm[0];
            double pmMemRemain = pm[1];

            double pmActualCpuUsed = pmActualResourceList.get(pmCount)[0];
            double pmActualMemUsed = pmActualResourceList.get(pmCount)[1];




            // If the remaining resource is enough for the container
            // And the OS is compatible
            if (pmCpuRemain >= vmCpuCapcacity &&
                    pmMemRemain >= vmMemCapacity) {

                Double pmScore = EvolveVmAllocationMethod(
                        input,
                        state,
                        ind,
                        threadnum,
                        stack,
                        pmCpuRemain,
                        pmMemRemain,
                        vmCpuCapcacity,
                        vmMemCapacity,
                        pmActualCpuUsed,
                        pmActualMemUsed,
                        vmUsedCpu,
                        vmUsedMem
                        );

                // Core of BestFit, score the bigger the better
                if (chosenPM == null || pmScore > BestScore) {
                    chosenPM = pmCount;
                    BestScore = pmScore;
                }
            } // End if
            // If there is no suitable PM (no PM has enough resources), then we just return null.
            pmCount ++;
        }
        return chosenPM;
    }

    private Double EvolveVmAllocationMethod(
                    final DoubleData input,
                    final EvolutionState state,
                    final Individual ind,
                    final int threadnum,
                    final ADFStack stack,
                    double pmCpuRemain,
                    double pmMemRemain,
                    double vmCpuCapacity,
                    double vmMemCapacity,
                    double pmActualCpuUsed,
                    double pmActualMemUsed,
                    double vmActualCpuUsed,
                    double vmActualMemUsed
    ){

        currentPmCpuRemain = pmCpuRemain;
        currentPmMemRemain = pmMemRemain;
        normalizedPmCpuRemain = currentPmCpuRemain / maxCPU;
        normalizedPmMemRemain = currentPmMemRemain / maxMem;
        normalizedVmCpuCapacity = vmCpuCapacity / maxCPU;
        normalizedVmMemCapacity = vmMemCapacity / maxMem;
        normalizedPmActualCpuUsed = pmActualCpuUsed / maxCPU;
        normalizedPmActualMemUsed = pmActualMemUsed / maxMem;
        normalizedVmActualCpuUsed = vmActualCpuUsed / maxCPU;
        normalizedVmActualMemUsed = vmActualMemUsed / maxMem;

        // update state in myEvolutionState
        MyEvolutionState myEvolutionState = (MyEvolutionState) state;
        myEvolutionState.currentPmCpuRemain = currentPmCpuRemain;
        myEvolutionState.currentPmMemRemain = currentPmMemRemain;
        myEvolutionState.normalizedPmCpuRemain = normalizedPmCpuRemain;
        myEvolutionState.normalizedPmMemRemain = normalizedPmMemRemain;
        myEvolutionState.normalizedVmCpuCapacity = normalizedVmCpuCapacity;
        myEvolutionState.normalizedVmMemCapacity = normalizedVmMemCapacity;
        myEvolutionState.normalizedPmActualCpuUsed = normalizedPmActualCpuUsed;
        myEvolutionState.normalizedPmActualMemUsed = normalizedPmActualMemUsed;
        myEvolutionState.normalizedVmActualCpuUsed = normalizedVmActualCpuUsed;
        myEvolutionState.normalizedPmActualMemUsed = normalizedPmActualMemUsed;

        // Evaluate the GP rule
        ((GPIndividual) ind).trees[0].child.eval(
                state, threadnum, input, stack, (GPIndividual) ind, containerAllocationProblem);
        return input.x;

    }


    private Double EvolveSelectionCreationMethod(
            DoubleData input,
            final EvolutionState state,
            final Individual ind,
            final int threadnum,
            final ADFStack stack,
            double vmCpuRemain,
            double vmMemRemain,
            double vmCpuCapacity,
            double vmMemCapacity,
            ArrayList<Double[]> pmResourceList,
            ArrayList<Double[]> actualPmResourceList){
        MyEvolutionState myEvolutionState = (MyEvolutionState)state;

        // allocated flag indicates whether the existing PM can host a newly created VM
        // true means YES, otherwise we must create new PM to host the VM
        boolean allocated = false;

//        DoubleData input = (DoubleData) (this.input);

        // The resource is normalized by the PM's capacity.
        normalizedContainerCpu = containerCpu / maxCPU;
        normalizedContainerMem = containerMem / maxMem;

        // update the data in myState

        myEvolutionState.normalizedVmCpuRemain = vmCpuRemain / maxCPU;
        myEvolutionState.normalizedVmMemRemain = vmMemRemain / maxMem;
        myEvolutionState.normalizedContainerCpu = normalizedContainerCpu;
        myEvolutionState.normalizedContainerMem = normalizedContainerMem;


        // we only consider the overhead of new VM
        if(newVmFlag) {
            normalizedVmCpuOverhead = vmCpuCapacity * vmCpuOverheadRate / maxCPU;
            normalizedVmMemOverhead = vmMemOverhead / maxMem;
        } else {
            normalizedVmCpuOverhead = 0;
            normalizedVmMemOverhead = 0;
        }

        // update the data in myState
        myEvolutionState.normalizedVmCpuOverhead = normalizedVmCpuOverhead;
        myEvolutionState.normalizedVmMemOverhead = normalizedVmMemOverhead;


        // Evaluate the GP rule
        ((GPIndividual) ind).trees[0].child.eval(
                state, threadnum, input, stack, (GPIndividual) ind, containerAllocationProblem);


        return input.x;
    }


    /**
     *
     * @return
     */
    private Integer VMSelectionCreation(
                                        final DoubleData input,
                                        final EvolutionState state,
                                        final Individual ind,
                                        final int threadnum,
                                        final ADFStack stack,
                                        ArrayList<Double[]> vmResourceList,
                                        ArrayList<Double[]> pmResourceList,
                                        ArrayList<Double[]> actualPmResourceList,
                                        ArrayList<Double[]> vmTypeList,
                                        HashMap<Integer, Integer> vmIndexTypeMapping,
                                        Double containerCpu,
                                        Double containerMem,
                                        int containerOS,
                                        double globalCpuWaste,
                                        double globalMemWaste
    ){
        Integer chosenVM = null;
        Double BestScore = null;
        int vmNum = vmResourceList.size();
        int vmCount = 0;
        ArrayList<Double[]> tempVMResourceList = new ArrayList<>();
        // make a copy of vmResourceList
        for (Double[] vmR: vmResourceList
             ) {
            tempVMResourceList.add(vmR.clone());
        }
//        tempVMResourceList = (ArrayList<Double[]>) vmResourceList.clone();
//        System.out.println(tempVMResourceList);
        for(Double[] vm:vmTypeList){
            // add this new VM into the tempList
            tempVMResourceList.add(new Double[]{
                    vm[0] - vm[0] * vmCpuOverheadRate,
                    vm[1] - vmMemOverhead,
                    new Double(containerOS)
            });
        }


        // Loop through the tempResourceList
        for(Double[] vm:tempVMResourceList){

            // Check if the vm exists
            newVmFlag = vmCount >= vmNum;

            // Get the remaining VM resources and OS
            double vmCpuRemain = vm[0];
            double vmMemRemain = vm[1];
            int vmOS = vm[2].intValue();
            int vmType ;
            if(vmCount < vmNum){
                vmType = vmIndexTypeMapping.get(vmCount);}
            else
                vmType = vmCount - vmNum;
                vmOS = containerOS;
            // If the remaining resource is enough for the container
            // And the OS is compatible
            if (vmCpuRemain >= containerCpu &&
                    vmMemRemain >= containerMem &&
                    vmOS == containerOS) {

                Double vmScore = EvolveSelectionCreationMethod(
                        input,
                        state,
                        ind,
                        threadnum,
                        stack,
                        vmCpuRemain,
                        vmMemRemain,
                        vmTypeList.get(vmType)[0],
                        vmTypeList.get(vmType)[1],
//                        globalCpuWaste,globalMemWaste,
                        pmResourceList,
                        actualPmResourceList);

                // Core of BestFit, score the bigger the better
                if (chosenVM == null || vmScore > BestScore) {
                    chosenVM = vmCount;
                    BestScore = vmScore;
                }

            } // End if

            // Increment the VM counter
            vmCount += 1;
        }

        newVmFlag = false;
        if (chosenVM ==null) {
            return vmResourceList.size()+1;}
        return chosenVM;

    }

}
//------------------------------------------------------------------PM creation with Ind---------------------------------------------------------------------------------
//package main;
//
//import com.opencsv.CSVReader;
//import ec.EvolutionState;
//import ec.Individual;
//import ec.gp.ADFStack;
//import ec.gp.GPIndividual;
//import ec.simple.SimpleEvolutionState;
////import sun.java2d.pipe.SpanShapeRenderer;
//
//import java.io.IOException;
//import java.io.Reader;
//import java.nio.file.Files;
//import java.nio.file.Paths;
//import java.sql.SQLOutput;
//import java.util.ArrayList;
//import java.util.HashMap;
//
//public class MainAllocationProcessContainer {
//    public double maxCPU = 211200;
//    public double maxMem = 4096000;
//    double maxVMMem = 128000;
//    public final double PMENERGY;
//
//
//
//    private final double vmCpuOverheadRate;
//    private final double vmMemOverhead;
//
//    public double containerCpu;
//    public double containerMem;
//    public int containerOs;
//    public double containerOsPro;
//
//
//    public double normalizedContainerCpu;
//    public double normalizedContainerMem;
//    public double normalizedVmCpuCapacity;
//    public double normalizedVmMemCapacity;
//
//
//    public double currentPmCpuRemain;
//    public double currentPmMemRemain;
//
//    public double normalizedVmCpuOverhead;
//    public double normalizedVmMemOverhead;
//
//    public double normalizedGlobalCpuWaste;
//    public double normalizedGlobalMemWaste;
//
//    public double normalizedPmCpuRemain;
//    public double normalizedPmMemRemain;
//
//    public double normalizedPmActualCpuUsed;
//    public double normalizedPmActualMemUsed;
//
//    public double normalizedVmActualCpuUsed;
//    public double normalizedVmActualMemUsed;
//
//
//
//    public boolean newVmFlag = false;
//
//
//    private final double k;
//
//
//    private ContainerAllocationProblem containerAllocationProblem;
//
//    public ArrayList<Integer> rentedPM = new ArrayList<>();
//    private boolean newPmFlag;
//
//    // Constructor
//    public MainAllocationProcessContainer(
//            ContainerAllocationProblem containerAllocationProblem,
//            EvolutionState state,
//            double PMENERGY,
//            double vmCpuOverheadRate,
//            double vmMemOverhead,
//            double k)
//    {
//
//        this.containerAllocationProblem = containerAllocationProblem;
//        this.PMENERGY = PMENERGY;
//        this.vmCpuOverheadRate = vmCpuOverheadRate;
//        this.vmMemOverhead = vmMemOverhead;
//        this.k = k;
//        this.rentedPM = new ArrayList<>();
//
//        // Initialize state
//
//
////        MyEvolutionState myEvolutionState = new MyEvolutionState();
//
//        MyEvolutionState myEvolutionState = (MyEvolutionState) state;
//
////        myEvolutionState.1.0 = 1.0;
////        myEvolutionState.1.0 = 1.0;
//        myEvolutionState.PMENERGY = PMENERGY;
//
//    }
//
//
//    private void initializeDataCenter(int testCase,
//                                      ArrayList<Double[]> pmResourceList,
//                                      ArrayList<Double[]> pmActualUsageList,
//                                      ArrayList<Double[]> vmResourceList,
//                                      ArrayList<Double[]> vmTypeList,
//                                      ArrayList<Double[]> pmTypeList,
//                                      ArrayList<ArrayList> initPm,
//                                      ArrayList<ArrayList> initVm,
//                                      ArrayList<ArrayList> initContainer,
//                                      ArrayList<ArrayList> initOs,
//                                      ArrayList<ArrayList> initPmType,
//                                      HashMap<Integer, Integer> VMPMMapping,
//                                      HashMap<Integer, Integer> vmIndexTypeMapping){
//
//        ArrayList<Double[]> initPmList = initPm.get(testCase);
//        ArrayList<Double[]> initVmList = initVm.get(testCase);
//        ArrayList<Double[]> containerList = initContainer.get(testCase);
//        ArrayList<Double[]> osList = initOs.get(testCase);
//        ArrayList<Double[]> initPmTypeList = initPmType.get(testCase);
//
//        int globalVmCounter = 0;
//        for(int i =0; i < initPmTypeList.size(); ++i){
//            int typePM =(initPmTypeList.get(i)[0]).intValue();
//            Double[] vms = initPmList.get(i);
//            double pmCPU = pmTypeList.get(typePM)[0];
//            double pmMem = pmTypeList.get(typePM)[1];
//            // Add the PM to resource List at the beginning stage
//            // The order of elements in the list is CPU, Memory and the type of this PM
//            pmResourceList.add(new Double[]{pmCPU,pmMem,(Double)(typePM*1.0)});
//            pmActualUsageList.add(new Double[]{pmCPU,pmMem});
//            // for this vm
//            for (int vmCounter = 0; vmCounter <vms.length; ++vmCounter ){
//                // Get the type of this VM
//                int vmType = vms[vmCounter].intValue() ;
//
//                // Get the OS type
//                Double[] os = osList.get(vmCounter + globalVmCounter);
//
//                // Create this VM
//                vmResourceList.add(new Double[]{
//                        vmTypeList.get(vmType)[0] - vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
//                        vmTypeList.get(vmType)[1] - vmMemOverhead,
//                        new Double(os[0])
//                });
//                // get the containers allocated on this VM
//                Double[] containers = initVmList.get(vmCounter + globalVmCounter);
//
//                // Allocate the VM to this PM,
//                // Allocation includes two part, first, pmResourceList indicates the left resource of PM (subtract entire VMs' size)
//                // pmIndex denotes the last PM. pmIndex should be at least 0.
//                int pmIndex = pmResourceList.size() - 1;
//
//                // update the pm left resources
//                pmResourceList.set(pmIndex, new Double[]{
//                        pmResourceList.get(pmIndex)[0] - vmTypeList.get(vmType)[0],
//                        pmResourceList.get(pmIndex)[1] - vmTypeList.get(vmType)[1],
//                        pmIndex*1.0
//                });
//
//                // The second part of allocation,
//                // We update the actual usage of PM's resources
//                pmActualUsageList.set(pmIndex, new Double[]{
//                        pmActualUsageList.get(pmIndex)[0] - vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
//                        pmActualUsageList.get(pmIndex)[1] - vmMemOverhead
//                });
//
//                // Map the VM to the PM
//                VMPMMapping.put(vmCounter + globalVmCounter, pmIndex);
//
//                // for each container
//                for(int conContainer = containers[0].intValue();
//                    conContainer < containers[containers.length - 1].intValue();
//                    ++conContainer){
//
//                    // Get the container's cpu and memory
//                    Double[] cpuMem = containerList.get(conContainer);
//
//                    //Create this container
//                    // get the left resources of this VM
//                    int vmIndex = vmResourceList.size() - 1;
//                    Double[] vmCpuMem = vmResourceList.get(vmIndex);
//
//                    // update the vm
//                    vmResourceList.set(vmIndex, new Double[] {
//                            vmCpuMem[0] - cpuMem[0],
//                            vmCpuMem[1] - cpuMem[1],
//                            new Double(os[0])
//                    });
//                    // Whenever we create a new VM, map its index in the VMResourceList to its type for future purpose
//                    for (int j = 0; j < vmResourceList.size(); j++) {
//                        vmIndexTypeMapping.put(j, vmResourceList.get(j)[2].intValue());
//                    }
////                    vmIndexTypeMapping.put(vmResourceList.size() - 1, vmType);
//                    // Add the Actual usage to the PM
//                    // Here, we must consider the overhead
//                    Double[] pmCpuMem = pmActualUsageList.get(pmIndex);
//
//                    // update the pm
//                    pmActualUsageList.set(pmIndex, new Double[]{
//                            pmCpuMem[0] - cpuMem[0],
//                            pmCpuMem[1] - cpuMem[1]
//                    });
//
//
////                    vm.addContainer(container);
//                } // Finish allocate containers to VMs
//            } // End  of each VM
//            // we must update the globalVmCounter
//            globalVmCounter += vms.length;
//
//        }// End  of each PM
//
////        System.out.println(globalVmCounter);
//    }
//
//
//
//
//    /**
//     * Calculate the energy consumption using the following equation:
//     * Energy = k * MaxEnergy + (1 - k)  * MaxEnergy * utilization_of_a_PM
//     * @param pmActualUsageList
//     * @return the energy consumption
//     */
//    private Double energyCalculation(ArrayList<Double[]> pmActualUsageList, ArrayList<Double[]> pmResourceList){
//        Double energy = 0.0;
//        for (int i = 0; i < pmResourceList.size(); i++) {
//            energy+= (pmActualUsageList.get(i)[0])/maxCPU *PMENERGY *(1-k) + k*PMENERGY * maxCPU;
//        }
////        for(Double[] pmActualResource:pmActualUsageList){
////            energy += ((PMCPU - pmActualResource[0]) / 1.0) * PMENERGY * (1 - k) + k * PMENERGY * this.maxCPU;
////        energy += ((PMCPU - pmActualResource[0]) / PMCPU) * PMENERGY * (1 - k) + k * PMENERGY;
////        }
//        return energy;
//    }
//
//
//    public ArrayList<Double> evaluate(
//            DoubleData input,
//            EvolutionState state,
//            ArrayList<ArrayList<Double[]>> inputX,
//            ArrayList<ArrayList> initVm,
//            ArrayList<ArrayList> initContainer,
//            ArrayList<ArrayList> initOs,
//            ArrayList<ArrayList> initPm,
//            ArrayList<ArrayList> initPmType,
//            ArrayList<Double[]> pmTypeList,
//            ArrayList<Double[]> vmTypeList,
//            GPIndividual vmSelectionCreationRule,
//            GPIndividual vmAllocationRule,
//            int threadnum,
//            ADFStack stack){
//
//        MyEvolutionState myEvolutionState = (MyEvolutionState) state;
//        // testCaseNum equals the current generation
//        int testCase = state.generation;
//        // initialize the resource lists
//        ArrayList<Double> resultList = new ArrayList<>();
//        ArrayList<Double> comparedResultList = new ArrayList<>();
//
//        // Loop through the testCases
////        for (int testCase = 0; testCase <= end - start - 1; ++testCase) {
//
//
//        double globalCPUWaste = 0;
//        double globalMEMWaste = 0;
//        ArrayList<Double[]> pmResourceList = new ArrayList<>();
//        ArrayList<Double[]> pmActualUsageList = new ArrayList<>();
//        ArrayList<Double[]> vmResourceList = new ArrayList<>();
//        HashMap<Integer, Integer> VMPMMapping = new HashMap<>();
//        HashMap<Integer, Integer> vmIndexTypeMapping = new HashMap<>();
//
//        // Initialize data center
//        initializeDataCenter(
//                testCase,
//                pmResourceList,
//                pmActualUsageList,
//                vmResourceList,
//                vmTypeList,
//                pmTypeList,
//                initPm,
//                initVm,
//                initContainer,
//                initOs,
//                initPmType,
//                VMPMMapping,
//                vmIndexTypeMapping);
//
//        // the total energy
//        Double Energy = energyCalculation(pmActualUsageList, pmResourceList);
////            Double pmFit = pmUtilRemain(pmActualUsageList, pmResourceList);
//
//
//        // No data center initialization
////                Double Energy = 0.0;
//        // Get all the containers
//        ArrayList<Double[]> containers = inputX.get(testCase);
//
//
//        int containerCounter = 0;
//        // Start simulation
//        for (Double[] container:containers) {
//
//            containerCpu = container[0];
//            containerMem = container[1];
//            containerOs = container[2].intValue();
//
//            // update myEvolutionState
//            myEvolutionState.containerCpu = containerCpu;
//            myEvolutionState.containerMem = containerMem;
//            myEvolutionState.containerOs = containerOs;
//
//
//
//
//            Integer chosenVM;
//            Integer currentVmNum = vmResourceList.size();
//            for (int i = 0; i < vmResourceList.size(); i++) {
//                vmIndexTypeMapping.put(i, vmResourceList.get(i)[2].intValue());
//            }
//
//            // select or create a VM
//            chosenVM = VMSelectionCreation(
//                    input,
//                    state,
//                    vmSelectionCreationRule,
//                    threadnum,
//                    stack,
//                    vmResourceList,
//                    pmResourceList,
//                    pmActualUsageList,
//                    vmTypeList,
//                    vmIndexTypeMapping,
//                    containerCpu,
//                    containerMem,
//                    containerOs,
//                    globalCPUWaste,
//                    globalMEMWaste);
//
//            // check if the VM exists, if chosenVM < currentVmNum is true, it means
//            // the chosenVM exists, we just need to update its resources
//
//            if (chosenVM < currentVmNum) {
//                // update the VM resources, allocating this container into this VM
//                vmResourceList.set(chosenVM, new Double[]{
//                        vmResourceList.get(chosenVM)[0] - containerCpu,
//                        vmResourceList.get(chosenVM)[1] - containerMem,
//                        new Double(containerOs)
//                });
//
//                // Find the pmIndex in the mapping
//                int pmIndex = VMPMMapping.get(chosenVM);
//
//                // update the PM actual resources
//                pmActualUsageList.set(pmIndex, new Double[]{
//                        pmActualUsageList.get(pmIndex)[0] - containerCpu,
//                        pmActualUsageList.get(pmIndex)[1] - containerMem
//                });
//
//                // Else, we need to create this new VM
//            } else {
//
//                // Retrieve the type of select VM
//                int vmType = chosenVM - currentVmNum;
//
//                // create this new VM
//                vmResourceList.add(new Double[]{
//                        vmTypeList.get(vmType)[0] - containerCpu - vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
//                        vmTypeList.get(vmType)[1] - containerMem - vmMemOverhead,
//                        new Double(containerOs)
//                });
//
//                // Whenever we create a new VM, map its index in the VMResourceList to its type for future purpose
//                vmIndexTypeMapping = new HashMap<>();
//                for (int i = 0; i < vmResourceList.size() ; i++) {
//                    vmIndexTypeMapping.put(i,vmResourceList.get(i)[2].intValue());
//                }
////                    vmIndexTypeMapping.put(vmResourceList.size() - 1, vmType);
//
//                // After creating a VM, we will choose a PM to allocate
//                Integer chosenPM = VMAllocation(
//                        input,
//                        myEvolutionState,
//                        vmAllocationRule,
//                        threadnum,
//                        stack,
//                        pmResourceList,
//                        pmActualUsageList,
//                        pmTypeList,
//                        vmTypeList.get(vmType)[0],
//                        vmTypeList.get(vmType)[1],
//                        containerCpu + vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
//                        containerMem + vmMemOverhead
//                );
//
//                // If we cannot choose a PM
//                if (chosenPM >= pmResourceList.size()) {
//
//                    // Add the VM to the newly created PM
//                    // We don't need to consider the overhead here.
//                    pmResourceList.add(new Double[]{
//                            pmTypeList.get(chosenPM-pmResourceList.size())[0] - vmTypeList.get(vmType)[0],
//                            pmTypeList.get(chosenPM-pmResourceList.size())[1] - vmTypeList.get(vmType)[1],
//                            chosenPM-pmResourceList.size()*1.0
//                    });
//
//                    // Add the Actual usage to the PM
//                    // Here, we must consider the overhead
//                    pmActualUsageList.add(new Double[]{
//                            pmTypeList.get(chosenPM-pmResourceList.size()+1)[0]- containerCpu - vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
//                            pmTypeList.get(chosenPM-pmResourceList.size()+1)[0] - containerMem - vmMemOverhead
//
//                    });
//
//                    // Map the VM to the PM
//                    VMPMMapping.put(vmResourceList.size() - 1, pmResourceList.size() - 1);
//
//
//                    // If there is an existing PM, we allocate it to an existing PM
//                } else {
//
//                    currentPmCpuRemain = pmResourceList.get(chosenPM)[0] - vmTypeList.get(vmType)[0];
//                    currentPmMemRemain = pmResourceList.get(chosenPM)[1] - vmTypeList.get(vmType)[1];
//                    double type =  pmResourceList.get(chosenPM)[2];
//                    // update the PM resources
//                    // pm resources - vm size
//                    pmResourceList.set(chosenPM, new Double[]{
//                            currentPmCpuRemain,
//                            currentPmMemRemain,
//                            type
//                    });
//
//                    // update the actual resources
//                    // Actual usage - container required - vm overhead
//                    pmActualUsageList.set(chosenPM, new Double[]{
//                            pmActualUsageList.get(chosenPM)[0] - containerCpu - vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
//                            pmActualUsageList.get(chosenPM)[1] - containerMem - vmMemOverhead
//                    });
//
//                    // Map the VM to the PM
//                    VMPMMapping.put(vmResourceList.size() - 1, chosenPM);
//
//                } // End of allocating a VM to an existing PM
//
//            } // End of creating a new VM
//
//
//            Energy = energyCalculation(pmActualUsageList, pmResourceList);
//
//            resultList.add(Energy);
//        } // End of all test cases
//
//
//        return resultList;
//    }
//    private Double[] pmCreation(ArrayList<Double[]> pmTypeList, double vmCpu, double vmMem) {
//        int chosedType = -1;
//        double bestFitValue = 0;
//        double bestcurrentUtil_CPU =0;
//        double bestcurrentUtil_Mem = 0;
//        double requireCPU = vmCpu;
//        double requireMem = vmMem;
//
//        for (int i = 0; i < pmTypeList.size(); i++) {
//            if(requireCPU < pmTypeList.get(i)[0] && requireMem < pmTypeList.get(i)[1]) {
//                double currentUtil_CPU = (pmTypeList.get(i)[0] -requireCPU ) / pmTypeList.get(i)[0];
//                double currentUtil_Mem = (pmTypeList.get(i)[1] - requireMem) / pmTypeList.get(i)[1];
//                double fitValue =  bestcurrentUtil_CPU * bestcurrentUtil_Mem;
//                if (bestcurrentUtil_CPU < currentUtil_CPU) {
//                    chosedType = i;
//                    bestcurrentUtil_CPU = currentUtil_CPU;
////                    bestcurrentUtil_Mem = currentUtil_Mem;
////                    bestFitValue = fitValue;
//                }
////                else if (fitValue > bestFitValue){// do not have both requirements better than others but together is better
////                    bestcurrentUtil_Mem  = currentUtil_Mem;
////                    bestcurrentUtil_CPU = currentUtil_CPU;
////                    bestFitValue = fitValue;
////                    chosedType = i;
////                }
//                else{//the current type is worse
//                    continue;
//                }
//            }
//
//        }
//        if (chosedType <0 || chosedType>pmTypeList.size()) chosedType =0;
//        Double newPM[] = new Double[]{pmTypeList.get(chosedType)[0], pmTypeList.get(chosedType)[1], chosedType*1.0};
//        return newPM;
//    }
//
//
//
////    /**
////     * VMAllocation allocates a VM to a PM
////     * @param pmResourceList is a list of remaining resources of PM, not the acutal used resources
////     * @param vmCpuCapacity the cpu capacity of a VM
////     * @param vmMemCapacity the memory capacity of a VM
////     * @return an index of a PM
////     */
////    private Integer VMAllocation(ArrayList<Double[]> pmResourceList, double vmCpuCapacity, double vmMemCapacity){
////        Integer chosenPM = null;
//////        Double BestScore = null;
////
////        // Loop through the PMs in the existing PM list
////        for(int pmCount = 0; pmCount < pmResourceList.size(); ++pmCount){
////            double pmCpuRemain = pmResourceList.get(pmCount)[0];
////            double pmMemRemain = pmResourceList.get(pmCount)[1];
////
////            // First Fit
////            if (pmCpuRemain >= vmCpuCapacity && pmMemRemain >= vmMemCapacity) {
////                chosenPM = pmCount;
////                break;
////            } // End if
////        }
////
////        return chosenPM;
////    }
//
//
//    /*
//     * Basically, VMAllocation still use the BestFit framework
//     */
//    private Integer VMAllocation(
//            DoubleData input,
//            final EvolutionState state,
//            final Individual ind,
//            final int threadnum,
//            final ADFStack stack,
//            ArrayList<Double[]> pmResourceList,
//            ArrayList<Double[]> pmActualResourceList,
//            ArrayList<Double[]> pmTypeList,
//            double vmCpuCapacity,
//            double vmMemCapacity,
//            double vmUsedCpu,
//            double vmUsedMem
//    ){
//
//        Integer chosenPM = null;
//        Double BestScore = null;
//        int pmCount = 0;
//        int pmNum = pmResourceList.size();
//        ArrayList<Double[]> tempPMResourceList = new ArrayList<>();
//        for (Double[] t : pmResourceList){
//            tempPMResourceList.add(t.clone());
//        }
//        // add the empty PMs
//        for(Double[] t_empty : pmTypeList){
//            tempPMResourceList.add(t_empty.clone());
//        }
//
//        // Loop through the tempResourceList
//        for(Double[] pm:tempPMResourceList){
//            newPmFlag = pmCount >= pmNum;
//            // Get the remaining PM resources
//            double pmCpuRemain = pm[0];
//            double pmMemRemain = pm[1];
//            double pmActualCpuUsed = 0;
//            double pmActualMemUsed = 0;
//            if(pmCount < pmNum){
//                pmActualCpuUsed = pmActualResourceList.get(pmCount)[0];
//                pmActualMemUsed = pmActualResourceList.get(pmCount)[1];
//            }
//            else{
//                pmActualCpuUsed = pm[0];
//                pmActualMemUsed = pm[1];
//            }
//
//            // If the remaining resource is enough for the container
//            // And the OS is compatible
//            if (pmCpuRemain >= vmCpuCapacity &&
//                    pmMemRemain >= vmMemCapacity) {
//
//                Double pmScore = EvolveVmAllocationMethod(
//                        input,
//                        state,
//                        ind,
//                        threadnum,
//                        stack,
//                        pmCpuRemain,
//                        pmMemRemain,
//                        vmCpuCapacity,
//                        vmMemCapacity,
//                        pmActualCpuUsed,
//                        pmActualMemUsed,
//                        vmUsedCpu,
//                        vmUsedMem
//                );
//
//                // Core of BestFit, score the bigger the better
//                if (chosenPM == null || pmScore > BestScore) {
//                    chosenPM = pmCount;
//                    BestScore = pmScore;
//                }
//            } // End if
//            // If there is no suitable PM (no PM has enough resources), then we just return null.
//            pmCount ++;
//        }
//        return chosenPM;
//    }
//
//    private Double EvolveVmAllocationMethod(
//            final DoubleData input,
//            final EvolutionState state,
//            final Individual ind,
//            final int threadnum,
//            final ADFStack stack,
//            double pmCpuRemain,
//            double pmMemRemain,
//            double vmCpuCapacity,
//            double vmMemCapacity,
//            double pmActualCpuUsed,
//            double pmActualMemUsed,
//            double vmActualCpuUsed,
//            double vmActualMemUsed
//    ){
//
//        currentPmCpuRemain = pmCpuRemain;
//        currentPmMemRemain = pmMemRemain;
//        normalizedPmCpuRemain = currentPmCpuRemain / maxCPU;
//        normalizedPmMemRemain = currentPmMemRemain / maxMem;
//        normalizedVmCpuCapacity = vmCpuCapacity / maxCPU;
//        normalizedVmMemCapacity = vmMemCapacity / maxMem;
//        normalizedPmActualCpuUsed = pmActualCpuUsed / maxCPU;
//        normalizedPmActualMemUsed = pmActualMemUsed / maxMem;
//        normalizedVmActualCpuUsed = vmActualCpuUsed / maxCPU;
//        normalizedVmActualMemUsed = vmActualMemUsed / maxMem;
//
//        // update state in myEvolutionState
//        MyEvolutionState myEvolutionState = (MyEvolutionState) state;
//        myEvolutionState.currentPmCpuRemain = currentPmCpuRemain;
//        myEvolutionState.currentPmMemRemain = currentPmMemRemain;
//        myEvolutionState.normalizedPmCpuRemain = normalizedPmCpuRemain;
//        myEvolutionState.normalizedPmMemRemain = normalizedPmMemRemain;
//        myEvolutionState.normalizedVmCpuCapacity = normalizedVmCpuCapacity;
//        myEvolutionState.normalizedVmMemCapacity = normalizedVmMemCapacity;
//        myEvolutionState.normalizedPmActualCpuUsed = normalizedPmActualCpuUsed;
//        myEvolutionState.normalizedPmActualMemUsed = normalizedPmActualMemUsed;
//        myEvolutionState.normalizedVmActualCpuUsed = normalizedVmActualCpuUsed;
//        myEvolutionState.normalizedPmActualMemUsed = normalizedPmActualMemUsed;
//
//        // Evaluate the GP rule
//        ((GPIndividual) ind).trees[0].child.eval(
//                state, threadnum, input, stack, (GPIndividual) ind, containerAllocationProblem);
//        return input.x;
//
//    }
//
//
//    private Double EvolveSelectionCreationMethod(
//            DoubleData input,
//            final EvolutionState state,
//            final Individual ind,
//            final int threadnum,
//            final ADFStack stack,
//            double vmCpuRemain,
//            double vmMemRemain,
//            double vmCpuCapacity,
//            double vmMemCapacity,
//            ArrayList<Double[]> pmResourceList,
//            ArrayList<Double[]> actualPmResourceList){
//        MyEvolutionState myEvolutionState = (MyEvolutionState)state;
//
//        // allocated flag indicates whether the existing PM can host a newly created VM
//        // true means YES, otherwise we must create new PM to host the VM
//        boolean allocated = false;
//
////        DoubleData input = (DoubleData) (this.input);
//
//        // The resource is normalized by the PM's capacity.
//        normalizedContainerCpu = containerCpu / maxCPU;
//        normalizedContainerMem = containerMem / maxMem;
//
//        // update the data in myState
//
//        myEvolutionState.normalizedVmCpuRemain = vmCpuRemain / maxCPU;
//        myEvolutionState.normalizedVmMemRemain = vmMemRemain / maxMem;
//        myEvolutionState.normalizedContainerCpu = normalizedContainerCpu;
//        myEvolutionState.normalizedContainerMem = normalizedContainerMem;
//
//
//        // we only consider the overhead of new VM
//        if(newVmFlag) {
//            normalizedVmCpuOverhead = vmCpuCapacity * vmCpuOverheadRate / maxCPU;
//            normalizedVmMemOverhead = vmMemOverhead / maxMem;
//        } else {
//            normalizedVmCpuOverhead = 0;
//            normalizedVmMemOverhead = 0;
//        }
//
//        // update the data in myState
//        myEvolutionState.normalizedVmCpuOverhead = normalizedVmCpuOverhead;
//        myEvolutionState.normalizedVmMemOverhead = normalizedVmMemOverhead;
//
//
//        // Evaluate the GP rule
//        ((GPIndividual) ind).trees[0].child.eval(
//                state, threadnum, input, stack, (GPIndividual) ind, containerAllocationProblem);
//
//
//        return input.x;
//    }
//
//
//    /**
//     *
//     * @return
//     */
//    private Integer VMSelectionCreation(
//            final DoubleData input,
//            final EvolutionState state,
//            final Individual ind,
//            final int threadnum,
//            final ADFStack stack,
//            ArrayList<Double[]> vmResourceList,
//            ArrayList<Double[]> pmResourceList,
//            ArrayList<Double[]> actualPmResourceList,
//            ArrayList<Double[]> vmTypeList,
//            HashMap<Integer, Integer> vmIndexTypeMapping,
//            Double containerCpu,
//            Double containerMem,
//            int containerOS,
//            double globalCpuWaste,
//            double globalMemWaste
//    ){
//        Integer chosenVM = null;
//        Double BestScore = null;
//        int vmNum = vmResourceList.size();
//        int vmCount = 0;
//        ArrayList<Double[]> tempVMResourceList = new ArrayList<>();
//        // make a copy of vmResourceList
//        for (Double[] vmR: vmResourceList
//        ) {
//            tempVMResourceList.add(vmR.clone());
//        }
////        tempVMResourceList = (ArrayList<Double[]>) vmResourceList.clone();
////        System.out.println(tempVMResourceList);
//        for(Double[] vm:vmTypeList){
//            // add this new VM into the tempList
//            tempVMResourceList.add(new Double[]{
//                    vm[0] - vm[0] * vmCpuOverheadRate,
//                    vm[1] - vmMemOverhead,
//                    new Double(containerOS)
//            });
//        }
//
//
//        // Loop through the tempResourceList
//        for(Double[] vm:tempVMResourceList){
//
//            // Check if the vm exists
//            newVmFlag = vmCount >= vmNum;
//
//            // Get the remaining VM resources and OS
//            double vmCpuRemain = vm[0];
//            double vmMemRemain = vm[1];
//            int vmOS = vm[2].intValue();
//            int vmType ;
//            if(vmCount < vmNum){
//                vmType = vmIndexTypeMapping.get(vmCount);}
//            else
//                vmType = vmCount - vmNum;
//            vmOS = containerOS;
//            // If the remaining resource is enough for the container
//            // And the OS is compatible
//            if (vmCpuRemain >= containerCpu &&
//                    vmMemRemain >= containerMem &&
//                    vmOS == containerOS) {
//
//                Double vmScore = EvolveSelectionCreationMethod(
//                        input,
//                        state,
//                        ind,
//                        threadnum,
//                        stack,
//                        vmCpuRemain,
//                        vmMemRemain,
//                        vmTypeList.get(vmType)[0],
//                        vmTypeList.get(vmType)[1],
////                        globalCpuWaste,globalMemWaste,
//                        pmResourceList,
//                        actualPmResourceList);
//
//                // Core of BestFit, score the bigger the better
//                if (chosenVM == null || vmScore > BestScore) {
//                    chosenVM = vmCount;
//                    BestScore = vmScore;
//                }
//
//            } // End if
//
//            // Increment the VM counter
//            vmCount += 1;
//        }
//
//        newVmFlag = false;
//        if (chosenVM ==null) {
//            return vmResourceList.size()+1;}
//        return chosenVM;
//
//    }
//
//}
