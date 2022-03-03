package main;

import java.sql.SQLOutput;
import java.util.ArrayList;
import java.util.HashMap;

public class initializationDataCenterForAll {
    private double currentEnergyUnitTime = 0.0;
    private double energyConsumption = 0.0;
    private double vmCpuOverheadRate = 0.1;



    private double vmMemOverhead = 200;
    private ArrayList<Double> maxCPUs = new ArrayList<>();
    private ArrayList<Double[]> pmResourceLsit= new ArrayList<>();
    private ArrayList<Double[]> pmActualUsageList = new ArrayList<>();
    private ArrayList<Double[]> vmResourceList = new ArrayList<>();
    private HashMap<Integer, Integer> VMPMMapping = new HashMap<>();
    private HashMap<Integer, Integer> vmIndexTypeMapping = new HashMap<>();
    public initializationDataCenterForAll(int testCase,
                                          ArrayList<Double[]> pmResourceList,
                                          ArrayList<Double[]> pmActualUsageList,
                                          ArrayList<Double[]> vmResourceList,
                                          ArrayList<Double[]> vmTypeList,
                                          ArrayList<Double[]> pmTypeList,
                                          ArrayList<ArrayList> initPm,
                                          ArrayList<ArrayList> initVm,
                                          ArrayList<ArrayList> initContainer,
//                                      ArrayList<Double[]> initContainer,
                                          ArrayList<ArrayList> initOs,
                                          ArrayList<ArrayList> initPmType,
                                          HashMap<Integer, Integer> VMPMMapping,
                                          HashMap<Integer, Integer> vmIndexTypeMapping){
        this.energyConsumption = 0;
        this.currentEnergyUnitTime = 0.0;
        ArrayList<Double[]> initPmList = initPm.get(testCase);
        ArrayList<Double[]> initVmList = initVm.get(testCase);
        ArrayList<Double[]> containerList = initContainer.get(testCase);
        ArrayList<Double[]> osList = initOs.get(testCase);
        ArrayList<Double[]> initPmTypeList = initPmType.get(testCase);
//        System.out.println("Initialization for allocation process container on ");

        int globalVmCounter = 0;
        // for each pm (in each pm it contains the type of vm holding by this pm)
        for(int i =0; i < initPmTypeList.size(); ++i){
            int typePM =(initPmTypeList.get(i)[0]).intValue();
            Double[] vms = initPmList.get(i);
            double pmCPU = pmTypeList.get(typePM)[0];
            double pmMem = pmTypeList.get(typePM)[1];
            double pmIdlePow = pmTypeList.get(typePM)[2];
            double pmMaxPow = pmTypeList.get(typePM)[3];
            double pmCoreNum = pmTypeList.get(typePM)[4];
            this.maxCPUs.add(new Double(pmCPU));
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
                        (vmTypeList.get(vmType)[0] *( 1-  vmCpuOverheadRate)),
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

                // update the pm left resources(bounds)
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
                        pmActualUsageList.get(pmIndex)[0] - (vmTypeList.get(vmType)[0] * vmCpuOverheadRate),
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
//                for(int conContainer = containers[0].intValue();
//                    conContainer < containers[containers.length - 1].intValue();
//                    ++conContainer){
                for(int conIndex = 0; conIndex < containers.length; conIndex++){
                    int conContainer = containers[conIndex].intValue();
                    // Get the container's cpu and memory
                    Double[] cpuMem = containerList.get(conContainer);

                    //Create this container
                    // get the left resources of this VM
                    int vmIndex = vmResourceList.size() - 1;
//                    if (vmIndex ==2){
//                        System.out.println("stop");
//                    }
                    Double[] vmCpuMem = vmResourceList.get(vmIndex);
                    double containerCPU = cpuMem[0];
                    double containerMem = cpuMem[1];
                    double vmCPULeft = vmCpuMem[0];
                    double vmMemLeft = vmCpuMem[1];
                    double checker = vmCPULeft - containerCPU;
                    vmMemLeft = vmMemLeft - containerMem;
                    // update the vm
//                    if (checker < 0){
//                        System.out.println("Wrong");
//                    }
                    vmResourceList.set(vmIndex, new Double[] {
                            vmCpuMem[0] - cpuMem[0],
                            vmCpuMem[1] - cpuMem[1],
//                            vmCPULeft ,
//                            vmMemLeft,
                            new Double(os[0]),
                            vmCpuMem[3] // num of core required by this vm
                    });
                    // Whenever we create a new VM, map its index in the VMResourceList to its type for future purpose
                    for (int j = 0; j < vmResourceList.size(); j++) {
                        vmIndexTypeMapping.put(j, vmResourceList.get(j)[2].intValue());
                    }
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
//                    if (pmActualUsageList.get(pmIndex)[0] < 0 ){
////                        System.out.println("in initialization pm usage wrong");
//                    }
                } // Finish allocate containers to VMs
            } // End  of each VM
            // we must update the globalVmCounter
            globalVmCounter += vms.length;

        }// End  of each PM

//        System.out.println(globalVmCounter);

        //Update the unit time energy consumption
        double unitPowerConsumption = 0.0;
        for(int i =0; i < pmResourceList.size(); i++){
//            Double type_d= pmResourceList.get(i)[2];
            Double idlePow = pmResourceList.get(i)[2];
            Double maxPow = pmResourceList.get(i)[3];
//            int type = type_d.intValue();
//            unitPowerConsumption += pmTypeList.get(type)[2]; //add the idle power to it
//            double utilization = (pmTypeList.get(type)[0] - pmActualUsageList.get(i)[0])/pmTypeList.get(type)[0];
            unitPowerConsumption += idlePow;
            double utilization = 1 - (pmActualUsageList.get(i)[0] /pmResourceList.get(i)[5]);
            if (utilization < 0) utilization = 1;
//            unitPowerConsumption += (pmTypeList.get(type)[3] - pmTypeList.get(type)[2])* (2 * utilization - Math.pow(utilization,1.4));
            unitPowerConsumption += (maxPow - idlePow) * (2*utilization - Math.pow(utilization, 1.4));
        }
        this.currentEnergyUnitTime = unitPowerConsumption;
        this.pmActualUsageList = pmActualUsageList;
        this.pmResourceLsit  = pmResourceList;
        this.vmResourceList = vmResourceList;
        this.vmIndexTypeMapping = vmIndexTypeMapping;
        this.VMPMMapping = VMPMMapping;
    }
    public double getUnitPower(){
//        this.currentEnergyUnitTime = 0.0;
        return this.currentEnergyUnitTime;
    }


    public double getEnergyConsumption() {
        return energyConsumption;
    }

    public double getVmCpuOverheadRate() {
        return vmCpuOverheadRate;
    }

    public double getVmMemOverhead() {
        return vmMemOverhead;
    }

    public ArrayList<Double[]> getPmResourceLsit() {
//        this.pmResourceLsit.clear();
        return pmResourceLsit;
    }

    public ArrayList<Double[]> getPmActualUsageList() {
//        this.pmActualUsageList.clear();
        return pmActualUsageList;
    }

    public ArrayList<Double[]> getVmResourceList() {
//        this.vmResourceList.clear();
        return vmResourceList;
    }

    public HashMap<Integer, Integer> getVMPMMapping() {
//        this.VMPMMapping.clear();
        return VMPMMapping;
    }

    public HashMap<Integer, Integer> getVmIndexTypeMapping() {
//        this.vmIndexTypeMapping.clear();
        return vmIndexTypeMapping;
    }
    public ArrayList<Double> getMaxCPUs() {
        return maxCPUs;
    }

}
