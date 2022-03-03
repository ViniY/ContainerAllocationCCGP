package main;
import ec.EvolutionState;
import ec.Individual;
import ec.gp.ADFStack;
import ec.gp.GPIndividual;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *  This file contains the allocation of VM to PM
 *  We newly changed the initializer function,
 *  the new initializer can read the PM and VM configuration into the PM resource list,
 *  it will be used as the initial state for the current generation,
 *  during the current generation the coming containers will be allocated to created VMs and the newly created VMs
 *  The newly created VMs will be allocated to desirable PM, these PM can be selected from already occupied PMs or
 *  we rent new PM from the PM type list.
 * **/
public class MainAllocationProcessVM_Ind {
    public double maxCPU = 211200; // Here we use the maximum value for both cpu and mem for normalization
    public double maxMem = 4096000;
//    public final double PMENERGY;

    // A list of containers
    private ArrayList<ArrayList<Double[]>> inputX = new ArrayList<>();

    // A list of candidate VMs, each vm has an array which includes its CPU and Mem capacity
    private ArrayList<Double[]> vmTypeList = new ArrayList<>();

    // An array of OS probability
    private ArrayList<Double> OSPro = new ArrayList<>();

    // An array of benchmark accumulated energy
    private ArrayList<Double> benchmark = new ArrayList<>();

    public boolean newPmFlag = false;
    // Initialization data
    ArrayList <ArrayList> initPm;
    ArrayList <ArrayList> initVm;
    ArrayList <ArrayList> initOs;
    ArrayList <ArrayList> initContainer;
    /**
    * Here are some Variable tracking the memory and cpu usage also some requirements about the container
    *
    * */
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

    public double normalizedPmCpuRemain;
    public double normalizedPmMemRemain;

    public double normalizedPmActualCpuUsed;
    public double normalizedPmActualMemUsed;

    public double normalizedVmActualCpuUsed;
    public double normalizedVmActualMemUsed;



    public boolean newVmFlag = false;

    private final double k;//overhead factror

    private ContainerAllocationProblem containerAllocationProblem;

    private double energyConsumption;
    public double currentEnergyUnitTime = 0;
    public double currentTimestamp=0.0;
    public double previousTimestamp = 0.0;
    private Double[] newPMHolder;
    private boolean newPMFlag;


    // Constructor
    public MainAllocationProcessVM_Ind(
                                    ContainerAllocationProblem containerAllocationProblem,
                                    EvolutionState state,
                                    double vmCpuOverheadRate,
                                    double vmMemOverhead,
                                    double k)
    {
        this.containerAllocationProblem = containerAllocationProblem; //from here we can access the VM types and PM types
        this.vmCpuOverheadRate = vmCpuOverheadRate;
        this.vmMemOverhead = vmMemOverhead;
        this.k = k;
        this.energyConsumption =0.0;
        this.currentEnergyUnitTime = 0.0;

        // Initialize state
        MyEvolutionState myEvolutionState = (MyEvolutionState) state;
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
        //TODO change
//        testCase = 0;
        this.currentEnergyUnitTime = 0.0;
        this.energyConsumption = 0.0;
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
            double pmIdlePow = pmTypeList.get(typePM)[2];
            double pmMaxPow = pmTypeList.get(typePM)[3];
            double pmCoreNum = pmTypeList.get(typePM)[4];
            // Add the PM to resource List at the beginning stage
            // The order of elements in the list is CPU, Memory and the type of this PM
            pmResourceList.add(new Double[]{pmCPU,pmMem, pmIdlePow, pmMaxPow, pmCoreNum, pmCPU});
            pmActualUsageList.add(new Double[]{pmCPU,pmMem, pmIdlePow, pmMaxPow, pmCoreNum, pmCPU});
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
                        new Double(os[0]),
                        vmTypeList.get(vmType)[2] // num of core required for this vm
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
                        pmResourceList.get(pmIndex)[1] - vmTypeList.get(vmType)[1],
                        pmResourceList.get(pmIndex)[2],
                        pmResourceList.get(pmIndex)[3],
                        pmResourceList.get(pmIndex)[4],
                        pmResourceList.get(pmIndex)[5]

                });
                // update the total cpu usage and total usage

                // The second part of allocation,
                // We update the actual usage of PM's resources
                pmActualUsageList.set(pmIndex, new Double[]{
                        pmActualUsageList.get(pmIndex)[0] - vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
                        pmActualUsageList.get(pmIndex)[1] - vmMemOverhead,
                        pmActualUsageList.get(pmIndex)[2],
                        pmActualUsageList.get(pmIndex)[3],
                        pmActualUsageList.get(pmIndex)[4],
                        pmActualUsageList.get(pmIndex)[5]

                });
                // update the actual usage

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
                            new Double(os[0]),
                            vmCpuMem[3] // num of core required by this vm
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
                            pmCpuMem[1] - cpuMem[1],
                            pmCpuMem[2],
                            pmCpuMem[3],
                            pmCpuMem[4],
                            pmCpuMem[5]

                    });
                    // update the total usage


//                    vm.addContainer(container);
                } // Finish allocate containers to VMs
            } // End  of each VM
            // we must update the globalVmCounter
            globalVmCounter += vms.length;

        }// End  of each PM


        //Update the unit time energy consumption at time zero
        double unitPowerConsumption = 0.0;
        for(int i =0; i < pmResourceList.size(); i++){
            Double idlePow = pmResourceList.get(i)[2];
            Double maxPow = pmResourceList.get(i)[3];
            unitPowerConsumption += idlePow;
            double utilization = 1 - (pmActualUsageList.get(i)[0] /pmResourceList.get(i)[5]);
            if (utilization < 0) System.out.println("Wrong utlization at initial stage ");
            unitPowerConsumption += (maxPow - idlePow) * (2*utilization - Math.pow(utilization, 1.4));
        }
        this.currentEnergyUnitTime = unitPowerConsumption;
    }

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

            ArrayList<Double[]> pmResourceList = new ArrayList<>();
            ArrayList<Double[]> pmActualUsageList = new ArrayList<>();
            ArrayList<Double[]> vmResourceList = new ArrayList<>();
            HashMap<Integer, Integer> VMPMMapping = new HashMap<>();
            HashMap<Integer, Integer> vmIndexTypeMapping = new HashMap<>();

            // Initialize data center
            initializeDataCenter(testCase,
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

            // No data center initialization
//                Double Energy = 0.0;
            // Get all the containers
//            ArrayList<Double[]> containers = inputX.get(testCase);
            ArrayList<Double[]> containers = inputX.get(0);
        // recording the first timestamp

        this.currentTimestamp = containers.get(0)[3];

            // Start simulation
            for (Double[] container:containers) {

                containerCpu = container[0];
                containerMem = container[1];
                containerOs = container[2].intValue();

                // update myEvolutionState
                myEvolutionState.containerCpu = containerCpu;
                myEvolutionState.containerMem = containerMem;
                myEvolutionState.containerOs = containerOs;

//                if (firstContainer){
//                    currentTimestamp = container[3];
//                    firstContainer = false;
//                }
                // if the timestamp changes we need to update the current power that consumed
                if (this.currentTimestamp - container[3]!=0){
                    boolean updatingEnergy = energyUpdate(container[3], currentTimestamp);
                    // update the timestamp
                    this.currentTimestamp = container[3];
                }

                Integer chosenVM;
                Integer currentVmNum = vmResourceList.size();

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
                        containerOs);

                // check if the VM exists, if chosenVM < currentVmNum is true, it means
                // the chosenVM exists, we just need to update its resources
                if (chosenVM < currentVmNum) {
                    // update the VM resources, allocating this container into this VM
                    vmResourceList.set(chosenVM, new Double[]{
                            vmResourceList.get(chosenVM)[0] - containerCpu,
                            vmResourceList.get(chosenVM)[1] - containerMem,
                            new Double(containerOs),
                            vmResourceList.get(chosenVM)[3]
                    });

                    // Find the pmIndex in the mapping
                    int pmIndex = VMPMMapping.get(chosenVM);
                    // update the unit power
                    if (containerCpu !=0){
                        double newRemain = pmActualUsageList.get(pmIndex)[0] - containerCpu;
                        updateCurrentPow(pmIndex, pmActualUsageList, newRemain);
                    }
                    // update the PM actual resources
                    pmActualUsageList.set(pmIndex, new Double[]{
                            pmActualUsageList.get(pmIndex)[0] - containerCpu,
                            pmActualUsageList.get(pmIndex)[1] - containerMem,
                            pmActualUsageList.get(pmIndex)[2],
                            pmActualUsageList.get(pmIndex)[3],
                            pmActualUsageList.get(pmIndex)[4],
                            pmActualUsageList.get(pmIndex)[5],
                    });
                    // Else, we need to create this new VM
                } else {

                    // Retrieve the type of select VM
                    int vmType = chosenVM - currentVmNum;

                    // create this new VM
                    vmResourceList.add(new Double[]{
                            vmTypeList.get(vmType)[0] - containerCpu - vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
                            vmTypeList.get(vmType)[1] - containerMem - vmMemOverhead,
                            new Double(containerOs),
                            vmTypeList.get(vmType)[2]
                    });

                    // Whenever we create a new VM, map its index in the VMResourceList to its type for future purpose
                    for (int i = 0; i < vmResourceList.size(); i++) {
                        vmIndexTypeMapping.put(i, vmResourceList.get(i)[2].intValue());

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
                                  containerMem + vmMemOverhead,
                                             vmTypeList.get(vmType)[2] // the core requirement of the vm
                                            );

                    // If we cannot choose a PM
                    if (chosenPM >= pmResourceList.size()) {
                        // find a new pm to create by best fit heuristic
//                        int chosedPMByBestFit = pmCreation(pmTypeList, vmTypeList.get(vmType)[0], vmTypeList.get(vmType)[1]);
//                        chosenPM = chosedPMByBestFit + pmResourceList.size();

                        // Add the VM to the newly created PM
                        // We don't need to consider the overhead here.
                        if (this.newPMFlag){
                            pmResourceList.add(this.newPMHolder.clone());
                            this.newPMFlag = false;
                        }


                        pmActualUsageList.add(new Double[]{this.newPMHolder[0] - containerCpu - vmTypeList.get(vmType)[0]* vmCpuOverheadRate,
                                this.newPMHolder[1] - containerMem - vmMemOverhead,
                                this.newPMHolder[2] ,
                                this.newPMHolder[3],
                                this.newPMHolder[4],
                                this.newPMHolder[5]});

                        // when adding a new pm we need to update the pow changes
                        this.currentEnergyUnitTime += pmActualUsageList.get(pmActualUsageList.size()-1)[2]; // add idle pow first
//                        double util = (containerCpu - vmTypeList.get(vmType)[0] * vmCpuOverheadRate)/ this.newPMHolder[0];
                        double util =1 -  (pmActualUsageList.get(pmActualUsageList.size()-1)[0] / pmActualUsageList.get(pmActualUsageList.size()-1)[5]);
                        if(util < 0) System.out.println("wrong util");
                        this.currentEnergyUnitTime += (this.newPMHolder[3] -this.newPMHolder[2]) * (2*util - Math.pow(util,1.4));
                        // Map the VM to the PM
                        VMPMMapping.put(vmResourceList.size() - 1, pmResourceList.size() - 1);
                        this.newPMFlag = false;

                        // If there is an existing PM, we allocate it to an existing PM
                    } else {

                        currentPmCpuRemain = pmResourceList.get(chosenPM)[0] - vmTypeList.get(vmType)[0];
                        currentPmMemRemain = pmResourceList.get(chosenPM)[1] - vmTypeList.get(vmType)[1];
                        if (containerCpu !=0) {
                            double newRemainingCPU = pmActualUsageList.get(chosenPM)[0] - containerCpu - vmTypeList.get(vmType)[0] * vmCpuOverheadRate;
                            updateCurrentPow(chosenPM, pmActualUsageList, newRemainingCPU);
                        }
                        // update the PM resources
                        // pm resources - vm size
                        pmResourceList.set(chosenPM, new Double[]{
                                currentPmCpuRemain,
                                currentPmMemRemain,
                                pmResourceList.get(chosenPM)[2],
                                pmResourceList.get(chosenPM)[3],
                                pmResourceList.get(chosenPM)[4],
                                pmResourceList.get(chosenPM)[5]
                        });

                        // update the actual resources
                        // Actual usage - container required - vm overhead
                        pmActualUsageList.set(chosenPM, new Double[]{
                                pmActualUsageList.get(chosenPM)[0] - containerCpu - vmTypeList.get(vmType)[0] * vmCpuOverheadRate,
                                pmActualUsageList.get(chosenPM)[1] - containerMem - vmMemOverhead,
                                pmActualUsageList.get(chosenPM)[2],
                                pmActualUsageList.get(chosenPM)[3],
                                pmActualUsageList.get(chosenPM)[4],
                                pmActualUsageList.get(chosenPM)[5],
                        });
                        // we aer using a created pm for vm allocation we need to update the total usage only

                        // Map the VM to the PM
                        VMPMMapping.put(vmResourceList.size() - 1, chosenPM);

                    } // End of allocating a VM to an existing PM

                } // End of creating a new VM


                // After each allocation of container, we evaluate the energy in the data-center
                // And add up the energy consumption to calculate the area under the curve
//                double increment = energyCalculation(pmActualUsageList);
//                Energy += increment;
//                pmFit += pmUtilRemain(pmActualUsageList, pmResourceList);
//                System.out.println(pmFit);

//            } // End of one testCase
//                ind.printIndividualForHumans(state, 0);
            // calculate the average increase of allocating a container for a taskCase,
            // We want to minimize this value
//            comparedResultList.add(Energy);
//                System.out.println("energy - sub = " + (Energy - benchmark.get(testCase)));
//            resultList.add(Energy / containers.size());
//                resultList.add(Energy);
//                resultList.add(pmFit);
//                resultList.add(energyCalculation(pmActualUsageList));

        } // End of all test cases

//        Double pmFit = pmRemain(pmResourceList);
//        resultList.add(pmFit);

        resultList.add(this.energyConsumption);
        this.energyConsumption =0.0;
        this.currentEnergyUnitTime = 0.0;
        return resultList;
    }


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
            double vmUsedMem,
            double vmCore
    ){

        Integer chosenPM = null;
        Double BestScore = null;
        int pmCount = 0;
        ArrayList<Double[]> tempPMResourceList = new ArrayList<>();

        for (Double[] t : pmResourceList){
            tempPMResourceList.add(t.clone());
        }
        // add the empty PMs
        for(Double[] t_empty : pmTypeList){
            if(t_empty[2] > vmCore){
                tempPMResourceList.add(t_empty.clone()); // only add the pm with capable core number
            }
        }
        // Loop through the tempResourceList
        for(Double[] pm:tempPMResourceList){
            // Get the remaining PM resources
            double pmCpuRemain = pm[0];
            double pmMemRemain = pm[1];
            double pmCore = pm[4];
            double pmActualCpuUsed = 0;
            double pmActualMemUsed = 0;
            int pmNum = pmResourceList.size();

            if(pmCount < pmNum){
                pmActualCpuUsed = pmActualResourceList.get(pmCount)[0];
                pmActualMemUsed = pmActualResourceList.get(pmCount)[1];
            }
            else{
                pmActualCpuUsed = pm[0];
                pmActualMemUsed = pm[1];
            }



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
                        vmUsedMem,
                        pmCore
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
        if(chosenPM >= pmResourceList.size()){
            this.newPMFlag = true;
            this.newPMHolder = new Double[]{tempPMResourceList.get(chosenPM)[0], // cpu
                    tempPMResourceList.get(chosenPM)[1], // mem
                    tempPMResourceList.get(chosenPM)[2], //
                    tempPMResourceList.get(chosenPM)[3],
                    tempPMResourceList.get(chosenPM)[4],
                    tempPMResourceList.get(chosenPM)[0]
            };
        }
        return chosenPM;
    }

    private Integer pmCreation(ArrayList<Double[]> pmTypeList, double vmCpu, double vmMem) {
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
//        Double newPM[] = new Double[]{pmTypeList.get(chosedType)[0], pmTypeList.get(chosedType)[1], chosedType*1.0};
        return chosedType;
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
                    double vmActualMemUsed,
                    double pmCore
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
        myEvolutionState.coreNumber = pmCore;

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
                                        int containerOS
    ){
        Integer chosenVM = null;
        Double BestScore = null;
        int vmNum = vmResourceList.size();
        int vmCount = 0;

        // make a copy of vmResourceList
        ArrayList<Double[]> tempVMResourceList = (ArrayList<Double[]>) vmResourceList.clone();
        for(Double[] vm:vmTypeList){
            // add this new VM into the tempList
            tempVMResourceList.add(new Double[]{
                    vm[0] - vm[0] * vmCpuOverheadRate,
                    vm[1] - vmMemOverhead,
                    new Double(containerOS),
                    vm[2]
            });
        }


        // Loop through the tempResourceList
        for(Double[] vm:tempVMResourceList){
            // Check if the vm existsw
            newVmFlag = vmCount >= vmNum;

            // Get the remaining VM resources and OS
            double vmCpuRemain = vm[0];
            double vmMemRemain = vm[1];
            int vmOS = vm[2].intValue();
            int vmType;
            if(vmCount < vmNum)
                vmType = vmResourceList.get(vmCount)[2].intValue();
            else
                vmType = vmCount - vmNum;

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
    public boolean updateCurrentPow(Integer chosenPM, ArrayList<Double[]> actualUsageList, double newCpuRemain){
        double privousUtil = (actualUsageList.get(chosenPM)[5] -  actualUsageList.get(chosenPM)[0])/  actualUsageList.get(chosenPM)[5];
        double maxPow = actualUsageList.get(chosenPM)[3];
        double minPow = actualUsageList.get(chosenPM)[2];
        double newUtil = (actualUsageList.get(chosenPM)[5] -  newCpuRemain)/  actualUsageList.get(chosenPM)[5];
        this.currentEnergyUnitTime  += (maxPow - minPow) * (2*newUtil - Math.pow(newUtil, 1.4) - 2*privousUtil + Math.pow(privousUtil, 1.4));
        return true;
    }
    /***
     * Update the power at current timestamp
     */

    private boolean energyUpdate(double curerntTime, double previousTime){
        this.energyConsumption += (this.currentEnergyUnitTime) * (curerntTime - previousTime) / 1000 / 3600;
        if (this.energyConsumption < 0){
            System.out.println("now is less than one");
        }
//        System.out.println("After updating the energy cost and now is (VM): " + this.energyConsumption);
        return true;

    }

}

