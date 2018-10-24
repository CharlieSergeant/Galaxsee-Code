using UnityEngine;
using System.Collections;
using System.Threading;


public abstract class TimestepModelAbstract : MonoBehaviour {

	public float modelDT = 0.01f;
	public Thread modelThread;

 	bool threadRunning = true;
	bool stepFree = true;
	bool stepRunning = false;
	public float modelT = 0.0f;
	float timescale = 1.0f;

	bool threaded = true;
	bool fastrun = true;
	bool paused = false;

	static readonly object _locker = new object();

	~TimestepModelAbstract() {
		Thread.Sleep (0);
		modelThread.Abort ();
	}

	public bool GetThreaded()
	{
		return threaded;
	}

	// Use this for initialization
	public void ModelStart()
	{
		modelT = 0.0f;
		// create thread before finishing Start
		if (threaded)
		{
			modelThread = new Thread(this.threadedActions);
			modelThread.Start();
		}

	}

	public void threadedActions()
	{
		while (threadRunning)
		{
			Thread.Sleep (0);

			try
			{
				if (stepFree || fastrun)
				{
					stepRunning = true;
					takeStep(modelDT*timescale);
					modelT += modelDT;//am I doing this twice?
					stepRunning = false;
				}
				if (!fastrun)
					lock (_locker)
					{
						stepFree = false;
					}
			}
			//(ThreadAbortException ex) 
			catch
			{
				threadRunning = false;
			}
		}
	}



	void FixedUpdate()
	{
		if (threaded)
		{
			lock (_locker)
			{
				stepFree = true;
			}
		}
		else
		{
			stepRunning = true;
			takeStep(modelDT);
			modelT += modelDT;
			stepRunning = false;
		}
	}

	public void Pause(bool yesNo)
	{
		lock (_locker)
		{
			stepFree = false;
			// would a lock be more efficient here?
			while (stepRunning) { } // wait for step to finish to avoid race condition
			// grow array if needed
			paused = yesNo;
		}
	}

	public abstract void takeStep (float dt);
}
