using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;


public class NBodyScript : Integrator
{

    int n;
    public double[] x;
    public double[] mass;
    public double softening_factor = 1.0e-3;
    public double sRad = 1.0e-3;
    public double G = 6.673e-11;
    public double m = 1.0;
    public double radius = 1.0;
    int start = 0;
    string readPath;
    string writePath;
    



    //file read in test --failed so far--
    /*public void fileSelected(string pathname)
    {
        string line;
        StreamReader sr = new StreamReader(pathname);
        while (IOExtras.ReadLine2(sr, out line))
        {
            int[] xyz = IOExtras.IntArray(line);
            x[start * 6 + 0] = xyz[0];
            x[start * 6 + 1] = xyz[1];
            x[start * 6 + 2] = xyz[2];
            x[start * 6 + 3] = xyz[3];
            x[start * 6 + 4] = xyz[4];
            x[start * 6 + 5] = xyz[5];
            mass[start] = xyz[6];
            start++;
        }
    }*/

    public void AllocNBS(int n)
    {
        this.n = n;
        Init(6 * n);
        x = new double[6 * n];
        mass = new double[n];
    }

    //this reads a file
    public void ReadFile(string fileName)
    {
        try
        {
            readPath = Application.dataPath + "/Resources/"+fileName+".txt";
            StreamReader reader = new StreamReader(readPath);
            string line;
            while (IOExtras.ReadLine2(reader, out line))
            {
                float[] xyz = IOExtras.FloatArray(line);
                x[start * 6 + 0] = (double)xyz[0];
                x[start * 6 + 1] = (double)xyz[1];
                x[start * 6 + 2] = (double)xyz[2];
                x[start * 6 + 3] = (double)xyz[3];
                x[start * 6 + 4] = (double)xyz[4];
                x[start * 6 + 5] = (double)xyz[5];
                mass[start] = (double)xyz[6];
                start++;
            }
        }
        catch
        {
            Debug.Log("Illegal argument in input file");
        }
    }

    public void WriteFile(bool isVirialCheck, double virial, float time)
    {
        try
        {
            writePath = Application.dataPath + "/Resources/output.txt";
            StreamWriter write = new StreamWriter(writePath, true);

            string line = "";
            //write.WriteLine("Next Frame");
            if (isVirialCheck)
            {
                line += virial+" "+time;
                write.WriteLine(line);
                line = "";
            }
            else
            {
                for (int i = 0; i < x.Length; i++)
                {
                    //formats string from scientific notation to decimal    
                    line += string.Format("{0:F20}", x[i].ToString("0.####################")) + " ";
                    //prints pos x,y,z and vel x,y,z and current time
                    if (i % 6 == 5 && (i != 0))
                    {
                        line += time;
                        write.WriteLine(line);
                        line = "";
                    }

                }
            }
            write.Close();
        }
        catch
        {
            Debug.Log("Cant write to file");
        }
    }


    public double setG(double mass,double length,double time)
    {
        return 6.6740831e-11 * mass * time * time / (length * length * length);
    }

    // Use this for initialization
    public void InitNBS()
    {
        //Allows for same starting position 
        Random.InitState(0);

        for (int i = start; i < n; i++)
        {
            bool done = false;
            while (!done)
            {
                x[i * 6 + 0] = Random.Range(-(float)radius, (float)radius); //pos x
                x[i * 6 + 1] = Random.Range(-(float)radius, (float)radius); //pos y
                x[i * 6 + 2] = Random.Range(-(float)radius, (float)radius); //pos z 
                double r2 = x[i * 6 + 0] * x[i * 6 + 0] +
                    x[i * 6 + 1] * x[i * 6 + 1] +
                    x[i * 6 + 2] * x[i * 6 + 2];
                if (r2 < radius * radius)
                    done = true;
            }
            //Debug.Log("x " + x[i * 6 + 0] + "y " + x[i * 6 + 1] + "z " + x[i * 6 + 0]);
            mass[i] = m / n; //initialize masses
            x[i * 6 + 3] = 0.0; //vel x
            x[i * 6 + 4] = 0.0; // vel y
            x[i * 6 + 5] = 0.0; //vel z
        }
        for (int i = start; i < n; i++)
        {
            double r2 = x[i * 6 + 0] * x[i * 6 + 0] +
                x[i * 6 + 1] * x[i * 6 + 1] +
                x[i * 6 + 2] * x[i * 6 + 2];
            double r = Mathd.Sqrt(r2);
            double v = 0.75 * Mathd.Sqrt(G * m * r / (radius * radius * radius));
            double th = Mathd.Atan2(x[i * 6 + 2], x[i * 6 + 0]);
            x[i * 6 + 3] = v * Mathd.Sin(th);
            x[i * 6 + 5] = -v * Mathd.Cos(th);
        }
    }

    public void InitNBS2()
    {
        //Allows for same starting position 
        Random.InitState(0);
        for (int i = start; i < n; i++)
        {
            bool done = false;
            while (!done)
            {
                x[i * 6 + 0] = Random.Range(-(float)radius, (float)radius); //pos x
                x[i * 6 + 1] = Random.Range(-(float)radius, (float)radius); //pos y
                x[i * 6 + 2] = Random.Range(-(float)radius, (float)radius); //pos z 
                double r2 = x[i * 6 + 0] * x[i * 6 + 0] +
                    x[i * 6 + 1] * x[i * 6 + 1] +
                    x[i * 6 + 2] * x[i * 6 + 2];
                if (r2 < radius * radius)
                    done = true;
            }
            //Debug.Log("x " + x[i * 6 + 0] + "y " + x[i * 6 + 1] + "z " + x[i * 6 + 0]);
            //mass[0] = 20;
            mass[i] = m / n; //initialize masses
        }
        double pe = getPE(x, mass);
        double v = Mathd.Sqrt(-pe / m);


        //Call getPE, because were going to use it to initialize the velocities.
        for (int j = start; j < n; j++)
        {
            Vector3 v3 = ((float)v) * Random.onUnitSphere;
            x[j * 6 + 3] = v3.x;
            x[j * 6 + 4] = v3.y;
            x[j * 6 + 5] = v3.z;

        }
    }

    public double getPE(double[] x, double[] mass)
    {
        /*Calculate Potential Energy of the System. Sum of triangle numbers of the system time complexity.
* Use this equation -> -G*((m1*m2)/r), where: G is G, m1 and m2 are masses, and r is the radius between them.*/
        double potentialEnergySum = 0.0;
        for (int i = 0; i < n; i++)
        {
            //Grab the mass of this particle.
            double m1 = mass[i];

            //Grab all the positions for the current point.
            double xCurrent = x[i * 6 + 0];
            double yCurrent = x[i * 6 + 1];
            double zCurrent = x[i * 6 + 2];

            /*Iterate all previous points of the system, and calculate the potential energy of each one,
             * against the current point. Add all those interactions together for the total potential energy
             * of just this CURRENT particle*/
            double currentPointsPotentialEnegy = 0.0;
            for (int j = i + 1; j < n; j++)
            {
                //Grab all the positions for this previous point.
                double xPrevious = x[j * 6 + 0];
                double yPrevious = x[j * 6 + 1];
                double zPrevious = x[j * 6 + 2];

                //Store the distances between this previous point and the current point for all x,y,z.
                double xDistance = xCurrent - xPrevious;
                double yDistance = yCurrent - yPrevious;
                double zDistance = zCurrent - zPrevious;

                //Calculate the distance between this previous point and the current point.
                double r = System.Math.Sqrt(xDistance * xDistance + yDistance * yDistance + zDistance * zDistance);

                //Calculate the potential energy of this interaction, and add it to THIS current points total PE.
                double m2 = mass[j];
                double thisPotentialEnergy = (-G) * ((m1 * m2) / r);
                currentPointsPotentialEnegy += thisPotentialEnergy;
            }
            //Add this points potential energ to the total.
            potentialEnergySum += currentPointsPotentialEnegy;
        }
        return potentialEnergySum;

    }

    public double getKE(double[] x, double[] mass)
    {
        /*Calculate Kinetic Energy of the System. Just the sum of all KE: (1/2)mv^2*/
        double kineticEnergySum = 0.0;
        for (int i = 0; i < n; i++)
        {
            double vx = x[i * 6 + 3];
            double vy = x[i * 6 + 4];
            double vz = x[i * 6 + 5];
            double totalVelocity = System.Math.Sqrt(vx * vx + vy * vy + vz * vz);

            kineticEnergySum += (0.5) * (mass[i] * totalVelocity * totalVelocity);
        }
        return kineticEnergySum;
    }

    override public void RatesOfChange(double[] x, double[] xdot)
    {
        for (int i = 0; i < n; i++)
        {
            xdot[i * 6 + 0] = x[i * 6 + 3];
            xdot[i * 6 + 1] = x[i * 6 + 4];
            xdot[i * 6 + 2] = x[i * 6 + 5];
            xdot[i * 6 + 3] = 0.0;
            xdot[i * 6 + 4] = 0.0;
            xdot[i * 6 + 5] = 0.0;
        }
        for (int i = 0; i < n; i++)
        {
            double xi = x[i * 6 + 0];
            double yi = x[i * 6 + 1];
            double zi = x[i * 6 + 2];
            for (int j = i + 1; j < n; j++)
            {
                double xj = x[j * 6 + 0];
                double yj = x[j * 6 + 1];
                double zj = x[j * 6 + 2];
                Vector3d drv = new Vector3d(xi - xj, yi - yj, zi - zj);
                double dr2 = (drv.sqrMagnitude + softening_factor * softening_factor);
                double dr = Mathd.Sqrt(dr2);
                Vector3d accel = -G / dr2 * drv / dr;
                xdot[i * 6 + 3] += accel.x * mass[j];
                xdot[i * 6 + 4] += accel.y * mass[j];
                xdot[i * 6 + 5] += accel.z * mass[j];
                xdot[j * 6 + 3] -= accel.x * mass[i];
                xdot[j * 6 + 4] -= accel.y * mass[i];
                xdot[j * 6 + 5] -= accel.z * mass[i];
            }
        }
    }
}
