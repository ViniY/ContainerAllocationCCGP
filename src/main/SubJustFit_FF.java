package main;

import ec.EvolutionState;
import functions.VMMEMOverhead;
import org.bytedeco.javacpp.presets.opencv_core;

import java.sql.SQLOutput;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * SubJustFit_FF represents that the algorithm contains three rules,
 * 1. Best-Fit with sub rule for allocating containers to VMs
 * 2. JustFit for VM creation, justFit finds the smallest VM for the container to create
 * 3. First-Fit for VM allocation
 */
public class SubJustFit_FF {
    private final ArrayList<ArrayList> initPmType;

    private double vmCpuOverheadRate;
    private double vmMemOverhead;
    private ArrayList<ArrayList> initVm;
    private ArrayList<ArrayList> initPm;
    private ArrayList<ArrayList> initOs;
    private ArrayList<ArrayList> initContainer;
    private ArrayList<Double[]> vmTypeList;
    private ArrayList<Double[]> pmTypeList;

    private ArrayList<Double> maxCPUs = new ArrayList<>();

    public double currentHour = 0.0;
    public double currentTimestamp;
    public double currentPower;
    public double currentEnergyConsumption;
    public boolean firstContainer = false;
    private ArrayList<ArrayList<Double[]>> inputX;

    public SubJustFit_FF(
            double vmCpuOverheadRate,
            double vmMemOverhead,
            ArrayList<ArrayList<Double[]>> inputX,
            ArrayList<ArrayList> initVm,
            ArrayList<ArrayList> initContainer,
            ArrayList<ArrayList> initOs,
            ArrayList<ArrayList> initPm,
            ArrayList<ArrayList> initPmType,
            ArrayList<Double[]> pmTypeList,
            ArrayList<Double[]> vmTypeList
    ) {
        this.pmTypeList = pmTypeList;
        this.vmCpuOverheadRate = vmCpuOverheadRate;
        this.vmMemOverhead = vmMemOverhead;
        this.initVm = initVm;
        this.initContainer = initContainer;
        this.initOs = initOs;
        this.initPm = initPm;
        this.inputX = inputX;
        this.vmTypeList = vmTypeList;
        this.initPmType = initPmType;
//    =====================================
        this.currentEnergyConsumption = 0.0;
        this.currentPower = 0.0;
        this.firstContainer = true;
    }


    public double allocate(
            int testCase
    ) {
        this.maxCPUs = new ArrayList<>();
        this.currentEnergyConsumption = 0.0;
        this.currentPower =0.0;
        ArrayList<Double[]> containers = inputX.get(testCase);
        // pmList, the real utilization of PM
        ArrayList<Double[]> pmList = new ArrayList<>();

        // pmStatusList, the VM boundary of PM
        ArrayList<Double[]> pmStatusList = new ArrayList<>();
//    private ArrayList<Integer[]> vmPmIndexMapping;
        ArrayList<Double[]> vmList = new ArrayList<>();
        HashMap<Integer, Integer> VMPMMapping = new HashMap<>();
        HashMap<Integer, Integer> vmIndexTypeMapping = new HashMap<>();
        //using a universal initializer
        initializationDataCenterForAll initializing = new initializationDataCenterForAll(
                testCase,
                pmStatusList,
                pmList,
                vmList,
                vmTypeList,
                pmTypeList,
                initPm,
                initVm,
                initContainer,
                initOs,
                initPmType,
                VMPMMapping,
                vmIndexTypeMapping);
        this.currentPower = initializing.getUnitPower();
        vmList = (ArrayList<Double[]>) initializing.getVmResourceList().clone();
        pmStatusList = (ArrayList<Double[]>) initializing.getPmResourceLsit().clone();
        pmList = (ArrayList<Double[]>) initializing.getPmActualUsageList().clone();
        vmIndexTypeMapping = (HashMap<Integer, Integer>) initializing.getVmIndexTypeMapping().clone();
        VMPMMapping = (HashMap<Integer, Integer>) initializing.getVMPMMapping().clone();
        this.maxCPUs = (ArrayList<Double>) initializing.getMaxCPUs().clone();
        vmList.clear();
        pmStatusList.clear();
        pmList.clear();
        vmIndexTypeMapping.clear();
        VMPMMapping.clear();
        this.currentPower = 0.0;
        this.maxCPUs.clear();
//        System.out.println("After initialization in heuristics process the power is " + this.currentPower);

        this.firstContainer = true;
        for (Double[] container : containers) {
            if (this.firstContainer) {
                this.currentTimestamp = container[3];
                this.firstContainer = false;
            }
            double containerCpu = container[0];
            double containerMem = container[1];
            int containerOs = container[2].intValue();
            double timestamp = container[3];
            if (timestamp != this.currentTimestamp) {
                updateEnergyConsumption(Math.abs(timestamp - this.currentTimestamp));
                this.currentTimestamp = timestamp; // update the current timestamp
            }
            if (!bestFit(container, pmList, vmList, VMPMMapping)) {
                int newVmIndex = justFit(container);
                // vm cpu, mem, os, core
                Double[] newVM = new Double[]{vmTypeList.get(newVmIndex)[0], vmTypeList.get(newVmIndex)[1], container[2], vmTypeList.get(newVmIndex)[2]};
                newVM[0] -= (containerCpu + newVM[0] * vmCpuOverheadRate); // deal the overhead for vm
                newVM[1] -= (containerMem + vmMemOverhead);
                vmList.add(newVM);
                vmIndexTypeMapping.put(vmList.size() - 1, newVmIndex);
//                System.out.println(vmList.size());
                if (!firstFit(newVmIndex, container, pmList, pmStatusList, vmList, VMPMMapping)) {
                    Double[] newPM = pmCreation(pmTypeList, vmTypeList.get(newVmIndex)[0], vmTypeList.get(newVmIndex)[1], vmTypeList.get(newVmIndex)[2]);
                    double maxCPU = newPM[0];
                    Double[] newPmStatus = newPM.clone();
                    newPM[0] -= (vmTypeList.get(newVmIndex)[0] * vmCpuOverheadRate );
                    newPM[0] -= containerCpu;
                    newPM[1] -= (vmMemOverhead + container[1]);
                    newPmStatus[0] -= vmTypeList.get(newVmIndex)[0];
                    newPmStatus[1] -= vmTypeList.get(newVmIndex)[1];
                    pmList.add(newPM);
                    this.maxCPUs.add(maxCPU);
                    pmStatusList.add(newPmStatus);

//                    vmPmIndexMapping.add(new Integer[]{vmList.size() - 1, pmList.size() - 1});
                    VMPMMapping.put(vmList.size() - 1, pmList.size() - 1);
//                    System.out.println(vmList.size() + " : " + pmList.size());
                    // Here we update the unit time power
                    updateCurrentPower(newPM, 0, maxCPU); //here we create a new pm
                    this.currentPower += newPM[2]; // add the idle power for the newly created pm
                }
            }
        }
        System.out.println("The final power unit time of heuristic is : " + this.currentPower);
        System.out.println("The power consumption is : " + this.currentEnergyConsumption);
        return this.currentEnergyConsumption;
    }



    private boolean firstFit(
            int vmIndex,
            Double[] container,
            ArrayList<Double[]> pmList,
            ArrayList<Double[]> pmStatusList,
            ArrayList<Double[]> vmList,
            HashMap<Integer, Integer> VMPMMapping) {
        boolean allocated = false;
        for (int i = 0; i < pmList.size(); i++) {
            Double[] pm = pmList.get(i);
            Double[] pmStatus = pmStatusList.get(i);
            double maxCPU = this.maxCPUs.get(i);
            double prevUtil = (maxCPU - pm[0]) / maxCPU;
            if (pmStatus[0] >= vmTypeList.get(vmIndex)[0] && // check cpu
                    pmStatus[1] >= vmTypeList.get(vmIndex)[1] && // check mem
                    pmStatus[4] >= vmTypeList.get(vmIndex)[2] && // check core
                    pm[0] >= (vmTypeList.get(vmIndex)[0] * vmCpuOverheadRate + container[0]) && // check container + overhead
                    pm[1] >= (vmMemOverhead + container[1]) // check containerMem + overhead mem, due to some container maybe super big
            ) {
                allocated = true;
                pm[0] -= (vmTypeList.get(vmIndex)[0] * vmCpuOverheadRate + container[0]);
                pm[1] -= (vmMemOverhead + container[1]);
                if(pm[0] < 0 ){
                    System.out.println("pm exceed limitation in firstFit");
                }
                pmStatus[0] -= vmTypeList.get(vmIndex)[0]; //update the bounds for pm
                pmStatus[1] -= vmTypeList.get(vmIndex)[1];
                VMPMMapping.put(vmList.size() - 1, i);
                updateCurrentPower(pm, prevUtil, maxCPU);

                break;
            }
        }
        return allocated;
    }

    // select the smallest (memory) VM type
    private int justFit(Double[] container) {
        int selectIndex = 0;
        int vmCounter = 0;
        double bestValue = Double.MAX_VALUE;
        for (Double[] vmType : vmTypeList) {
            if ((vmType[0]*(1-vmCpuOverheadRate))>=container[0] && (vmType[1] >=(container[1]+vmMemOverhead))){
                double leftMem = vmType[1] - container[1] - vmMemOverhead;
                if (leftMem < bestValue) {
                    selectIndex = vmCounter;
                    bestValue = leftMem;
                }
            }
            vmCounter++;
        }
        return selectIndex;
    }

    private boolean bestFit(
            Double[] container,
            ArrayList<Double[]> pmList,
            ArrayList<Double[]> vmList,
            HashMap<Integer, Integer> VMPMMapping) {

        boolean allocated = false;
        int vmCounter = 0;
        int bestVmIndex = 0;
        double bestValue = Double.MAX_VALUE;
        for (Double[] vm : vmList) {
            double vmCpu = vm[0];
            double vmMem = vm[1];
            int vmOs = vm[2].intValue();

            double containerCpu = container[0];
            double containerMem = container[1];
            int containerOs = container[2].intValue();

            if (vmCpu >= containerCpu && vmMem >= containerMem && vmOs == containerOs) {
                allocated = true;
                double cpuLeft = vm[0] - container[0];
                double memLeft = vm[1] - container[1];
                double subValue = Math.abs(cpuLeft - memLeft);
                if (subValue < bestValue) {
                    bestValue = subValue;
                    bestVmIndex = vmCounter;
                }
            }
            vmCounter++;
        }

        if (allocated) {
            Double[] vm = vmList.get(bestVmIndex);
            vm[0] -= container[0];
            vm[1] -= container[1];

            int pmIndex = VMPMMapping.get(bestVmIndex);
            Double[] pm = pmList.get(pmIndex);
            double maxCPU = this.maxCPUs.get(pmIndex);
            double prevUtil = (maxCPU - pm[0]) / maxCPU;
            pm[0] -= container[0];
            pm[1] -= container[1];
            if(pm[0] < 0 ){
                System.out.println(" Warning PM exceed limitation");
            }
            updateCurrentPower(pm, prevUtil, maxCPU);
        }
        return allocated;
    }

    // For PM creation create pm
    private Double[] pmCreation(ArrayList<Double[]> pmTypeList, double vmCpu, double vmMem, double vmCore) {
        int chosedType = -1;
        double bestFitValue = 0;
        double bestcurrentUtil_CPU = 0;
        double bestcurrentUtil_Mem = 0;
        double requireCPU = vmCpu;
        double requireMem = vmMem;
        for (int i = 0; i < pmTypeList.size(); i++) {
            if (requireCPU <= pmTypeList.get(i)[0] &&
                    requireMem <= pmTypeList.get(i)[1] &&
                    vmCore <= pmTypeList.get(i)[4] ) {
                double currentUtil_CPU = (requireCPU) / pmTypeList.get(i)[0]; // used / maxCPU
                if (bestcurrentUtil_CPU < currentUtil_CPU) {
                    chosedType = i;
                    bestcurrentUtil_CPU = currentUtil_CPU;
                }
                else {//the current type is worse
                    continue;
                }
            }

        }
//       // The pm type holding cpu, mem, idle, max, core num
        Double newPM[] = new Double[]{pmTypeList.get(chosedType)[0], pmTypeList.get(chosedType)[1], pmTypeList.get(chosedType)[2],
                pmTypeList.get(chosedType)[3], pmTypeList.get(chosedType)[4]};
        return newPM;
    }

    /***
     * Update the current energy consumption
     * @param  timeDif (difference between the new timestamp and the previous one)
     */

    public void updateEnergyConsumption(double timeDif) {
        this.currentEnergyConsumption += timeDif / 1000 / 3600 * this.currentPower;
    }

    public void updateCurrentPower(Double[] updatedPM, double prevUtil, double maxCPU) {
//        this.currentPower += updatedPM[2];
        double currentUtil = (maxCPU - updatedPM[0]) / maxCPU;
        double powToAdd = 0.0;
        if (currentUtil < 0 || currentUtil > 1) {
            if (currentUtil < 0) {
                System.out.println("The current util is below zero");
            }
            System.out.println("util greater than one");
        }
            if (prevUtil > currentUtil) {
                System.out.println("current util is less than before");
            }
            if (prevUtil == 0) {
                this.currentPower += (updatedPM[3] - updatedPM[2]) * (2 * currentUtil - Math.pow(currentUtil, 1.4));
            } else {

                powToAdd = (updatedPM[3] - updatedPM[2]) * (2 * currentUtil - Math.pow(currentUtil, 1.4) -
                        (2 * prevUtil - Math.pow(prevUtil, 1.4)));
            }
//            if (powToAdd < 0) System.out.println("Pow Decrease ? subjustFF");
            this.currentPower += powToAdd;

    }
}
