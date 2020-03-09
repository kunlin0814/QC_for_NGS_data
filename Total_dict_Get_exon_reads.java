import java.io.IOException;
import java.awt.List;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Collections;   

// The java version of the target_reads.py scripts. Each sample takes around 15 min vs 2 hr in python script
public class Total_dict_Get_exon_reads {


	public static void main(String[] args) throws IOException {
		
		String ExonInterval = "/work/szlab/kh31516_Lab_Share_script/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list";
        String SamFile = args[0];
        String file_name = args[1];
        int read_length = Integer.parseInt(args[3]);
        File file_output = new File(args[2]);	
		
		//Map<String, ArrayList<ExonLocationInfo>> small_interval_dict = new TreeMap<>();
		//Map<String, ArrayList<ExonLocationInfo>> large_interval_dict = new TreeMap<>();
		Map<String, ArrayList<ExonLocationInfo>> Total_interval_dict = new HashMap<>();
		
		BufferedReader objReader = null;
		objReader = new BufferedReader(new FileReader(ExonInterval));
		String strCurrentLine;
		while ((strCurrentLine = objReader.readLine()) != null) {
        ArrayList <ExonLocationInfo> ExonLocation = new ArrayList<ExonLocationInfo>(); 
            
        String[] tokens = strCurrentLine.split(":");
        String chrom = tokens[0];
        int start = Integer.parseInt(tokens[1].split("-")[0]);
        int end = Integer.parseInt(tokens[1].split("-")[1]);
       
        ExonLocation.add(new ExonLocationInfo(start, end));
        	if (Total_interval_dict.containsKey(chrom)== false) {
        		
        		Total_interval_dict.put(chrom, ExonLocation);
        	}
        	else {
        		Total_interval_dict.get(chrom).add(new ExonLocationInfo(start, end));
        		
        	}
        		
    
                
        }
		objReader.close();
		Set <String> Total_chrom_list =  Total_interval_dict.keySet();
		
		for (  String ele : Total_chrom_list) {
			Collections.sort(Total_interval_dict.get(ele));
		}
		
		
	ArrayList<String> summary = getExonCount(SamFile,read_length, Total_interval_dict);
	
	BufferedWriter bw = null;
	
	FileWriter fw = new FileWriter(file_output);
	bw = new BufferedWriter(fw);
	//bw.write("File_name\t"+ "Total_reads\t"+"Total_uniq\t"+"uniq_mapped_rate\t"+"Total_read_pairs\t"+"uniq_Exonic_region\t"+"uniq_Exonic_region_paris_rates\t"+'\n');
	bw.write(file_name+'\t'+ summary.get(0)+'\t'+summary.get(1)+'\t'+summary.get(2)+'\t'+summary.get(3)+'\t'+summary.get(4)+'\t'+summary.get(5)+'\n');
	bw.close();
	}
	
	public static ArrayList<String>getExonCount(String SamFile, int read_length, Map<String, ArrayList<ExonLocationInfo>> Total_interval_dict) throws NumberFormatException, IOException{
		int unique = 0 ;// 3 is the remainder of 16 in the flag
		//duplicate = 0 #3
		//Onemapped = 0 #5,9
		//incorrect = 0 #1
		//unmapped = 0 #13   
	    int total = 0 ;
		int pass_line = 0 ;
	
		ArrayList <String> transcript_list = new ArrayList<String>(); 
		ArrayList <String> summary_list = new ArrayList<String>(); 
		BufferedReader objReader_samfile = null;
		objReader_samfile = new BufferedReader(new FileReader(SamFile));
		String strCurrentLine_sam;
		while ((strCurrentLine_sam = objReader_samfile.readLine()) != null) {
			
		if(strCurrentLine_sam.startsWith("@")) {
			 ;
		}
		else {
			String [] file_lst =  strCurrentLine_sam.split("\t");
			total += 1 ;
		    String reads_name = file_lst[0];
            String reads_chr = 	file_lst[2];
            int reads_position = Integer.parseInt(file_lst[3]);
            int status = Integer.parseInt(file_lst[1])%16;
            if (status == 3) {
		            
		 if (file_lst[11].contains("XT:")) { 
            String status2 = file_lst[11].split(":")[2];
        if (status2.equals("U") || status2.equals("M")) {
            unique += 1 ;
       
        ArrayList<ExonLocationInfo> Total_exome_loc = Total_interval_dict.get(reads_chr);
        int location_status = binarySearch(Total_exome_loc, 0, (Total_exome_loc.size()-1),reads_position,read_length);
        if (location_status != -2) {
          transcript_list.add(reads_name);
           pass_line+=1;
    }
            }        
		    }
	}
		}      
		    
		}
		
		
	double uniq_mapped_rate = Double.valueOf(unique)/Double.valueOf(total);
	double uniq_Exonic_region_mapped_rate = Double.valueOf(pass_line)/Double.valueOf(unique);
	
	int pairs = total/2;
	summary_list.add(Integer.toString(total));
	summary_list.add(Integer.toString(unique));
	summary_list.add(String.valueOf(uniq_mapped_rate));
	summary_list.add(Integer.toString(pairs));
	summary_list.add(Integer.toString(pass_line));
	summary_list.add(String.valueOf(uniq_Exonic_region_mapped_rate));
	objReader_samfile.close();	
	
		
		return summary_list;
		
	}
	
	
    public static int withinRegion(int reads_position, int exom_start, int exom_end, int read_length) {
		int reads_end = reads_position + read_length-1 ;
    	int reads_start = reads_position ; 
    	int value = 0;
    	
    	// means the read is greater
    	if (reads_start > exom_end) {
    		value = 0; }
    	// means the read is smaller
    	else if (reads_end < exom_start) {
    		value = -1 ;}
    	// means the read is in the region
    	else {
    		value = 1;
    	}
     return value;
	
    }
    
    public static int binarySearch (ArrayList<ExonLocationInfo> arr, int left, int right, int reads_position, int read_length) { 
    	
		// Check base case 
        if (right >= left) {
        	
        	int mid = (left + right) / 2; 
            int exon_start_value = arr.get(mid).getstart();
            int exon_end_value = arr.get(mid).getend();
            // If element is present at the middle itself 
            int region_status = withinRegion( reads_position, exon_start_value, exon_end_value,  read_length );
            if (region_status == 1) {
                return mid ;
            }
            // If element is smaller than mid, then it can only be present in left subarray 
            else if (region_status  == -1){
                return binarySearch(arr, left, mid-1, reads_position, read_length) ;
                }
      
            // Else the element can only be present in right subarray 
            else if (region_status == 0){
                return binarySearch(arr, mid+1, right, reads_position, read_length);
                } 
        }
        	
            // Element is not present in the array 
        	return -2; 
    
    }
}
	
class ExonLocationInfo implements Comparable <ExonLocationInfo>{

    int start;
    int end;
    
    public ExonLocationInfo( int start, int end) {
        this.start = start;
        this.end = end;
    }
    public int getstart() {         
        return start;     
      }       
      public int getend() {         
        return end;     
      }                      
  public int compareTo(ExonLocationInfo exon) {              
	  // use start point to sort
	  return this.start-exon.start;     
    	  }                



}

