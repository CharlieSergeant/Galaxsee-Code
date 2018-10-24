using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ConsoleLog : MonoBehaviour {

	Text Output = null;
	InputField Input = null;
	string outputString = "";
	string inputString = "";
	ArrayList inputLines = null;

	// Use this for initialization
	void Awake () {
		Output = transform.Find ("Output").Find ("Container").Find ("OutputText").gameObject.GetComponent<Text>();
		Input = transform.Find ("Input").Find ("InputField").gameObject.GetComponent<InputField>();
		outputString = Output.text;
		inputLines = new ArrayList ();
	}

	public void Write(string line) {
		outputString = outputString + line;
		Output.text = outputString;
	}

	public string getLine() {
		if (inputLines.Count > 0) {
			string line = inputLines [0] as string;
			inputLines.RemoveAt (0);
			return line;
		} else {
			return null;
		}
	}

	public void processLine(string line) {
		inputString = Input.text;

		outputString = outputString + "> " + line + inputString + "\n";
		Output.text = outputString;
		Input.text = "";
		inputLines.Add (inputString);
		Scrollbar sb = Output.transform.parent.parent.Find ("Scrollbar").gameObject.GetComponent<Scrollbar>() as Scrollbar;
		sb.value = 0.0f;
	}
	
	// Update is called once per frame
	void Update () {
		
	}
}
