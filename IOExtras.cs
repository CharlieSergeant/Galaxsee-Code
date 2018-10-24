using UnityEngine;
using System.Collections;
using System.IO;

public class IOExtras
{

    /**
	 * ReadLine2(StreamReader sr, out string line) is a whitespace stripping, 
	 * blank line skipping, comment ignoring version of readline
	 * 
	 * as input it takes a streamreader and a string output variable
	 * it returns true or false based on whether a line was there to read.
	 * 
	 **/

    // need to strip double spaces


    public static char[] commentChar = new char[] { '#' };
    public static char[] whitespaceChar = new char[] { ' ', '\t' };
    public static char[] delimChar = new char[] { ' ', '\t', ',', ';' };

    public static string stripDoubleSpaces(string line)
    {
        foreach (char c in delimChar)
        {
            string pattern = c + " ";
            while (line.IndexOf(pattern) > -1)
            {
                line = line.Replace(pattern, c.ToString());
            }
            pattern = " " + c;
            while (line.IndexOf(pattern) > -1)
            {
                line = line.Replace(pattern, c.ToString());
            }
        }

        // fix errors of space with delimit character
        return line;
    }

    public static bool ReadLine2(StreamReader sr, out string line)
    {
        if (sr.Peek() >= 0)
        {
            while (true)
            {
                string testLine = sr.ReadLine();
                testLine = testLine.Trim();
                string[] halves = testLine.Split(commentChar, 2);
                testLine = halves[0].Trim();
                testLine = stripDoubleSpaces(testLine);
                if (testLine.Equals(""))
                {
                    if (sr.Peek() < 0)
                    {
                        line = null;
                        return false;
                    }
                }
                else
                {
                    line = testLine;
                    return true;
                }
            }
        }
        else
        {
            line = null;
            return false;
        }
    }


    public static int[] IntArray(string line)
    {
        string[] words = line.Split(delimChar);
        int[] values = new int[words.Length];
        for (int i = 0; i < words.Length; i++)
        {
            values[i] = int.Parse(words[i]);
        }
        return values;
    }

    public static float[] FloatArray(string line)
    {
        string[] words = line.Split(delimChar);
        float[] values = new float[words.Length];
        for (int i = 0; i < words.Length; i++)
        {
            values[i] = float.Parse(words[i]);
        }
        return values;
    }

    public static string[] StringArray(string line)
    {
        string[] words = line.Split(delimChar);
        return words;
    }

    public static int KeyVal(string line, out string key, out string val)
    {
        string[] words = line.Split(delimChar, 2);
        int n = words.Length;
        if (n > 0)
            key = words[0];
        else
            key = null;
        if (n > 1)
            val = words[1];
        else
            val = null;
        return n;
    }
}
