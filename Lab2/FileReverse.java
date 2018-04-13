//-----------------------------------------------------------------------------
// FileReverse.java
// Illustrates file IO and tokenization of strings.
//-----------------------------------------------------------------------------
import java.io.*;
import java.util.Scanner;

class FileReverse{
   public static void main(String[] args) throws IOException{

      int lineNumber = 0;

      // check number of command line arguments is at least 2
      if(args.length < 2){
         System.out.println("Usage: FileCopy <input file> <output file>");
         System.exit(1);
      }

      // open files
      Scanner in = new Scanner(new File(args[0]));
      PrintWriter out = new PrintWriter(new FileWriter(args[1]));

      // read lines from in, extract and print tokens from each line
      while( in.hasNextLine() ){
         lineNumber++;

         // trim leading and trailing spaces, then add one trailing space so 
         // split works on blank lines
         String line = in.nextLine().trim() + " "; 
//         out.println(line + "11 ");

         // split line around white space 
         String[] token = line.split("\\s+");  

         // print out tokens       
         int n = token.length;
//         out.println("Line " + lineNumber + " contains " + n + " tokens:");
         for(int i=0; i<n; i++){
            out.print(stringReverse(token[i],token[i].length()));
            out.println("");
         }
      }

      // close files
      in.close();
      out.close();
   }
   
      public static String stringReverse(String a, int n){
      if(n>0){
         String b=(a.substring(n-1,n));
         return (b + (stringReverse(a, n-1)));
      }
      return ("");
        
   }
   
}
