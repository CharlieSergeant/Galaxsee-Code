using UnityEngine;
using System.Collections;
using System.Runtime.InteropServices;
using System;
using UnityEngine.UI;


public class Model : TimestepModelAbstract {
    //namespace for all of the DLL imports in a C library

	[DllImport("UnityDLLExample",EntryPoint="setMass")]
	private static extern void setMass (IntPtr p, double value);

	[DllImport("UnityDLLExample",EntryPoint="setMasses")]
	private static extern void setMasses (IntPtr p, double [] value);

	[DllImport("UnityDLLExample",EntryPoint= "setSofteningNBodyModel")]
	private static extern void setSoftening (IntPtr p, double value);

    /*Shield radius methods*/
    [DllImport("UnityDLLExample", EntryPoint = "setSradNbodyModel")]
    private static extern void setSrad(IntPtr p, double value);
    /**/

    [DllImport("UnityDLLExample",EntryPoint="setX")]
	private static extern void setX (IntPtr p, double value, int i);

	[DllImport("UnityDLLExample",EntryPoint="setXArray")]
	private static extern void setXArray (IntPtr p, double [] value);

    [DllImport("UnityDLLExample", EntryPoint = "setTStep")]
    private static extern void setTStep(IntPtr p, double dt);

    [DllImport("UnityDLLExample", EntryPoint = "setG")]
    private static extern void setG(IntPtr p, double g);

    [DllImport("UnityDLLExample",EntryPoint="getX")]
	private static extern double getX (IntPtr p,int i);

    [DllImport("UnityDLLExample",EntryPoint="getXArray")]
	private static extern IntPtr getXArray (IntPtr p);

	[DllImport("UnityDLLExample",EntryPoint="allocModel")]
	private static extern IntPtr allocModel (int i);

	[DllImport("UnityDLLExample",EntryPoint="initModel")]
	private static extern void initModel (IntPtr p);

	[DllImport("UnityDLLExample",EntryPoint="stepModelRK2")]
	private static extern void stepModelRK2 (IntPtr p, double dt);

	[DllImport("UnityDLLExample",EntryPoint="stepModelRK4")]
	private static extern void stepModelRK4 (IntPtr p, double dt);

	[DllImport("UnityDLLExample",EntryPoint="stepModelEuler")]
	private static extern void stepModelEuler (IntPtr p, double dt);

    //Allows for the change of G based on scale of project
    /*[DllImport("UnityDLLExample", EntryPoint = "getGFromSI")]
    private static extern double getGFromSI(double mass, double length, double dt);*/

    //Leapfrog integration test
    [DllImport("UnityDLLExample", EntryPoint = "stepNbodyModelLeapfrog")]
    private static extern void stepNbodyModelLeapfrog(IntPtr p, double dt);

    //ABM integration test
    [DllImport("UnityDLLExample", EntryPoint = "stepNbodyModelABM")]
    private static extern void stepNbodyModelABM(IntPtr p, double dt);

    NBodyScript theModel;

    //Setting values for G scale
    public double setGMass = 5.972 * Math.Pow(10, 24); //Earth mass
    public double setGLength = 1.49 * Math.Pow(10, 11);//AU (distance between earth and sun)
    public double setGTime = 3.154*Math.Pow(10,7);  //Seconds in year

	public int n=1000;
    public int integrationVal = 1;
    GameObject [] theObjects;
	GameObject theParent;
	IntPtr pluginModel;
	double [] marshalledX;
    

    //public string integerationMethod = "RK4";
    public bool useDLL = true;
    public bool WriteToFile = true;
    public bool drawPixels = true;
	public bool drawSpheres = true;
    public bool trailsOn = true;
    public bool isVirialCheck = false;

	public float sphereSize = 0.1f;
	public Material PCMat;
	public Material sphereMat;

    //add softening potentials here 
	public float softening_factor = 1.0e-3f;
    public float sRad = 1.0e-3f;
    public float totalMass = 1.0f;

    //UI references
    public Text objectsText;
    public Text timeStepText;
    public Text SoftFacText;
    public Text SoftRadText;
    public Text ReadFileText;
    public Text TotalMasstext;
    public Text StarSizeText;

    GameObject ResetMainCamera;
    Camera camera;
    //Vector3 startPostion;
    public string fileName = "tester";
    public int forceVal = 1;

    private void Start()
    {
        Debug.Log("Model should not have started");
        modelT = 0;
        //startPostion= transform.position;
        //camera = GameObject.Find("BasicPlayer").GetComponent<Camera>();
        //startPostition= camera.transform.position;
        
    }

    void ModelReset()
    {
        Debug.Log("Resetting Model");
        Start();
        Destroy(theParent);
        useDLL = true;
        theModel = new NBodyScript();
        theModel.AllocNBS(n);
        if ((fileName.Equals("")))
        {

        }
        else
        {
            theModel.ReadFile(fileName);
        }
        theModel.m = totalMass;
        Debug.Log("G constant used is " + theModel.setG(setGMass, setGLength, setGTime));
        theModel.G = theModel.setG(setGMass, setGLength, setGTime);
        theModel.softening_factor = softening_factor;
        theModel.sRad = sRad;

        theModel.InitNBS2(); //initial creation of positions and velocities of system
        if (useDLL)
        { // goes to C code 
            pluginModel = allocModel(n);  //create model and allocate memory in C
            marshalledX = new double[6 * n];
            setSoftening(pluginModel, theModel.softening_factor); //softening factor
            setSrad(pluginModel, theModel.sRad); //softening radius
            setG(pluginModel, theModel.G);
            //setG (pluginModel, getGFromSI(setGMass,setGSize,setGTime)); //Sets G based on scale
            setMasses(pluginModel, theModel.mass);
            setXArray(pluginModel, theModel.x);
            //initModel (pluginModel);
        }

        theObjects = new GameObject[n];
        theParent = new GameObject();
        for (int i = 0; i < n; i++)
        {
            theObjects[i] = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            theObjects[i].GetComponent<Collider>().enabled = false;
            theObjects[i].GetComponent<Renderer>().material = sphereMat;
            if (trailsOn)
            {
                TrailRenderer tr = theObjects[i].AddComponent<TrailRenderer>();
                tr.material = PCMat;
                tr.material.color = Color.white;
                tr.minVertexDistance = 0.1f;
                tr.startWidth = 0.025f;
                tr.time = 120;
            }
            theObjects[i].transform.parent = theParent.transform;
        }
        /*
        if (useDLL) {
            for (int i = 0; i < n; i++) {
                bool inside = false;
                while (!inside) {
                    float R = 2.0f;
                    marshalledX [i * 6 + 0] = (double)UnityEngine.Random.Range (-R, R);
                    marshalledX [i * 6 + 1] = (double)UnityEngine.Random.Range (-R, R);
                    marshalledX [i * 6 + 2] = (double)UnityEngine.Random.Range (-R, R);
                    double x = marshalledX [i * 6 + 0];
                    double y = marshalledX [i * 6 + 1];
                    double z = marshalledX [i * 6 + 2];
                    double r = Mathd.Sqrt (x * x + y * y + z * z);
                    if (r <= R)
                        inside = true;
                    double v = 1.0;
                    marshalledX [i * 6 + 3] = -v * z / r;
                    marshalledX [i * 6 + 4] = 0.0f;
                    marshalledX [i * 6 + 5] = v * x / r;
                }

            }


            setXArray (pluginModel, marshalledX, 6 * n);
        }
    */
        ModelStart();
    }

    //UI
    public void onMenuTrigger(string val)
    {
        //Start Model Menu
        if (val == "start")
        {
            ModelReset();
        }
        //Objects Menu
        if (val == "objects_custom")
        {
            Debug.Log(objectsText.text+ " is the N value passed");
            n = Int32.Parse(objectsText.text);
        }
        if(val == "objects_1000")
        {
            Debug.Log("1000 objects");
            n = 1000;
        }
        if (val == "objects_100")
        {
            Debug.Log("100 objects");
            n = 100;
        }
        if(val == "objects_50")
        {
            Debug.Log("50 objects");
            n = 50;
        }
        if (val == "objects_25")
        {
            Debug.Log("25 objects");
            n = 25;
        }
        if (val == "objects_10")
        {
            Debug.Log("10 objects");
            n = 10;
        }
        //Time Step Menu
        if (val == "timeStep_custom")
        {
            Debug.Log(timeStepText.text + " is the TS value passed");
            modelDT = float.Parse(timeStepText.text);
        }
        if (val == "timeStep_0.01")
        {
            Debug.Log("Time Step set to 0.01");
            modelDT = 0.01f;
        }
        if (val == "timeStep_0.001")
        {
            Debug.Log("Time Step set to 0.001");
            modelDT = 0.001f;
        }
        if (val == "timeStep_0.0001")
        {
            Debug.Log("Time Step set to 0.0001");
            modelDT = 0.0001f;
        }
        if (val == "timeStep_0.00001")
        {
            Debug.Log("Time Step set to 0.00001");
            modelDT = 0.00001f;
        }
        //Softening Factor Menu
        if (val == "MenuSoft_custom")
        {
            Debug.Log(timeStepText.text + " is the TS value passed");
            softening_factor = float.Parse(SoftFacText.text);
        }
        if (val == "MenuSoft_0.1")
        {
            Debug.Log("Softening Factor set to 0.1");
            softening_factor = 0.1f;
        }
        if (val == "MenuSoft_0.01")
        {
            Debug.Log("Softening Factor set to 0.01");
            softening_factor = 0.01f;
        }
        if (val == "MenuSoft_0.001")
        {
            Debug.Log("Softening Factor set to 0.001");
            softening_factor = 0.001f;
        }
        if (val == "MenuSoft_0.0001")
        {
            Debug.Log("Softening Factor set to 0.0001");
            softening_factor = 0.0001f;
        }
        //Softening Radius Menu
        if (val == "MenuSoftRad_custom")
        {
            Debug.Log(timeStepText.text + " is the sRad value passed");
            sRad = float.Parse(SoftRadText.text);
        }
        if (val == "MenuSoftRad_0.1")
        {
            Debug.Log("Softening Radius set to 0.1");
            sRad = 0.1f;
        }
        if (val == "MenuSoftRad_0.01")
        {
            Debug.Log("Softening Radius set to 0.01");
            sRad = 0.01f;
        }
        if (val == "MenuSoftRad_0.001")
        {
            Debug.Log("Softening Radius set to 0.001");
            sRad = 0.001f;
        }
        if (val == "MenuSoftRad_0.0001")
        {
            Debug.Log("Softening Radius set to 0.0001");
            sRad = 0.0001f;
        }
        //Integration Menu
        if(val == "Integrate_Euler")
        {
            Debug.Log("Itegrating with Euler's Method");
            integrationVal = 1;
        }
        if (val == "Integrate_RK2")
        {
            Debug.Log("Itegrating with RK2");
            integrationVal = 2;
        }
        if (val == "Integrate_RK4")
        {
            Debug.Log("Itegrating with RK4");
            integrationVal = 3;
        }
        if (val == "Integrate_Bashforth")
        {
            Debug.Log("Itegrating with ABM");
            integrationVal = 4;
        }
        if (val == "Integrate_Leapfrog")
        {
            Debug.Log("Itegrating with Leapfrog Method");
            integrationVal = 5;
        }
        //Force Calc Menu
        if(val == "Direct_Force")
        {
            forceVal = 1;
        }
        if (val == "PPPM")
        {
            forceVal = 2;
        }
        if(val == "Barnes_Hut")
        {
            forceVal = 3;
        }
        //Mass of System Menu
        if(val == "MassSys_100")
        {
            totalMass = 100;
        }
        if (val == "MassSys_1000")
        {
            totalMass = 1000;
        }
        if (val == "MassSys_10000")
        {
            totalMass = 10000;
        }
        if (val == "MassSys_50000")
        {
            totalMass = 50000;
        }
        if (val == "MassSys_custom")
        {
            Debug.Log(timeStepText.text + " is the sRad value passed");
            totalMass = float.Parse(TotalMasstext.text);
        }
        //sphere size menu
        if(val == "StarSize_0.1")
        {
            sphereSize = 0.1f;
        }
        if (val == "StarSize_0.2")
        {
            sphereSize = 0.1f;
        }
        if (val == "StarSize_0.5")
        {
            sphereSize = 0.1f;
        }
        if (val == "StarSize_1")
        {
            sphereSize = 0.1f;
        }
        if (val == "StarSize_custom")
        {
            Debug.Log(StarSizeText.text + " is the sRad value passed");
            sphereSize = float.Parse(StarSizeText.text);
        }
        //Tracer On/Off
        if(val == "Tracer_On")
        {
            trailsOn = true;
        }
        if(val == "Tracer_Off")
        {
            trailsOn = false;
        }
        //Virial On/Off
        if (val == "Virial_On")
        {
            isVirialCheck = true;
        }
        if (val == "Virial_Off")
        {
            isVirialCheck = false;
        }
        //GScale
        if(val == "EarthMasses_AU")
        {
            theModel.G = theModel.setG(setGMass, setGLength, setGTime);
        }
        //initial shape
        if(val == "rings")
        if(val == "sphere")
        if(val == "cube")
        if(val == "2Dgrid")
        if(val == "3Dgrid")
        if(val == "disk")
        
        //Play/Pause
        if(val =="Play/Pause")
        {

        }

        //ResetCamera
        if (val == "ResetCamera")
        {

        }
        //Read File Menu
        if (val == "File_custom")
        {
            Debug.Log(ReadFileText.text + " is the custom file value passed");
            fileName = timeStepText.text;
        }
        if (val == "File_SolarSystem")
        {
            Debug.Log("Using file SolarSystem");
            fileName = "SolarSystem";
        }
        if (val == "File_EarthandSun")
        {
            Debug.Log("Using file Earth and Sun");
            fileName = "EarthandSun";
        }
        if (val == "File_NoVelocity")
        {
            Debug.Log("Using file No Velocity");
            fileName = "NoVelocity";
        }
        if(val == "3")
        {

        }
        else 
        {

        }
    }

    // Update is called once per frame
    void Update ()
    {
        if (theModel != null)
        {
            double[] x = new double[6 * n];
            if (useDLL)
            {
                for (int i = 0; i < 6 * n; i++)
                {
                    x[i] = getX(pluginModel, i);
                }
                //IntPtr foou = getXArray(pluginModel);
                //Marshal.Copy (foou, marshalledX, 0, 6*n);
                //x = marshalledX;
            }
            else
            {
                x = theModel.x;
            }

            /*Virial Theorem: The Kinetic Energy should be (-1/2)*PotentialEnergy.
                * To test this, adding the Kinetic Energy to 1/2 the Potential Energy should equal 0*/
            double kineticEnergy = theModel.getKE(x, theModel.mass);
            double potentialEnergy = theModel.getPE(x, theModel.mass);
            double virial = kineticEnergy + (0.5) * potentialEnergy;
            //Debug.Log("Virial: " + virial);


            Vector3[] test = new Vector3[n];
            Color[] cols = new Color[n];
            for (int i = 0; i < n; i++)
            { //changes the game objects position
                theObjects[i].transform.position = new Vector3(
                    (float)x[i * 6 + 0],
                    (float)x[i * 6 + 1],
                    (float)x[i * 6 + 2]);
                test[i] = theObjects[i].transform.position;
                Vector3 v = new Vector3((float)x[i * 6 + 3], (float)x[i * 6 + 4], (float)x[i * 6 + 5]);

                //write file every frame
                /* double[] xLine = new double[6];
                    xLine[0] = x[i * 6 + 0];
                    xLine[1] = x[i * 6 + 1];
                    xLine[2] = x[i * 6 + 2];
                    xLine[3] = x[i * 6 + 3];
                    xLine[4] = x[i * 6 + 4];
                    xLine[5] = x[i * 6 + 5];
                    theModel.WriteFile(xLine);*/
                if (WriteToFile)
                    theModel.WriteFile(isVirialCheck, virial, modelT);

                //color scale stuff
                float cscale = v.magnitude;
                cols[i] = new Color(cscale, 1.0f, 1.0f - cscale);
                theObjects[i].GetComponent<Renderer>().material.color = cols[i];
                if (trailsOn)
                {

                    theObjects[i].GetComponent<TrailRenderer>().startColor = cols[i];
                    theObjects[i].GetComponent<TrailRenderer>().endColor = cols[i];
                }
                theObjects[i].transform.localScale = sphereSize * Vector3.one;

                /*if (trailsOn)
                {
                    theObjects[i].GetComponent<TrailRenderer>().startWidth = 0.005f;
                    theObjects[i].GetComponent<TrailRenderer>().endWidth = 0.005f;
                }*/

                if (!drawSpheres) theObjects[i].SetActive(false);
                else theObjects[i].SetActive(true);
            }
            if (drawPixels)
            {
                Mesh mesh = CreatePointMesh(test);
                mesh.colors = cols;
                Graphics.DrawMesh(mesh, Vector3.zero,
                    Quaternion.identity, PCMat, 0);
            }
        }
	}

	Mesh CreatePointMesh(Vector3[] points)
	{
        
            Mesh mesh = new Mesh();
            mesh.vertices = points;
            // You can also apply UVs or vertex colours here.

            int[] indices = new int[points.Length];
            for (int i = 0; i < points.Length; i++)
                indices[i] = i;

            mesh.SetIndices(indices, MeshTopology.Points, 0);


            return mesh;
	}

	override public void takeStep (float dt)
	{
        if (useDLL)
        {
            if(integrationVal == 1)
                stepModelEuler(pluginModel, (double)dt);
            if (integrationVal == 2)
                stepModelRK2(pluginModel, (double)dt);
            if (integrationVal == 3)
                stepModelRK4(pluginModel, (double)dt);
            if (integrationVal == 4)
                stepNbodyModelABM(pluginModel, (double)dt);
            if (integrationVal == 5)
                stepNbodyModelLeapfrog(pluginModel, (double)dt);
        }
        else
        {
            theModel.RK4Step(theModel.x, dt);

        }
        modelT += dt;
    }	
}

