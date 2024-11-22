#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "htslib/sam.h"
void write_wig(int tid,int* cover,bam_hdr_t *header,FILE* wf);
void write_cover_tsv(int tid,int* cover,bam_hdr_t *header,FILE* wf);
double sample_number(int tid,int* cover,bam_hdr_t*header);
double t_mean(int tid,int* cover,bam_hdr_t*header);
double t_sigma(int tid,int* cover,bam_hdr_t*header);
double t_skewness(int tid,int* cover,bam_hdr_t*header);
int main(int argc, char **argv){//need the directory of bamfile as input
    bam_hdr_t *header;//record the header of the bam file
    bam1_t *aln = bam_init1();//keep the information of every row of bam file
    struct transcript*tlist;//record transcript information
    FILE*wf,*wt,*wc;
    samFile *in = sam_open(argv[1], "r");  //open bam file
    wf=fopen(argv[2],"a+");//create wig file
    //wt=fopen(argv[3],"a+");//create tsv file
    wc=fopen(argv[4],"a+");//create tsv2 file
    //fprintf(wt,"ENST\tlength\tread_count(bp)\tmean\tstandard_deviation\tskewness\n");
    header = sam_hdr_read(in);
    int stat=sam_read1(in, header, aln);
    /*
        record the status of input
        >= 0    read succesfully
        = -1    EOF
        < -1    error
    */
    int ctid=aln->core.tid; //record current transcript id
    int*tcover;//record reads coverage of current transcript
    tcover =(int*)malloc(sizeof(int)*(*(header->target_len+ctid)));
    for(int i=0;i<(*(header->target_len+ctid));i++)//initialization
        *(tcover+i)=0;        
    int flag=1;//use to break loop
    if(stat<-1) flag=0;
    if(aln->core.isize>0){//count coverage
        for(int i=(aln->core.pos);i<(aln->core.pos+aln->core.isize);i++){
                *(tcover+i)+=1;
            }
    }
    while(flag==1){
        stat=sam_read1(in, header, aln);
        if(stat<-1) flag=0;
        if(stat==-1) flag=2;
        if((aln->core.tid)!=ctid){
            write_wig(ctid,tcover,header,wf);
            write_cover_tsv(ctid,tcover,header,wc);
            //fprintf(wt,"%s\t%d\t%f\t%f\t%f\t%f\n",*(header->target_name+ctid),*(header->target_len+ctid),sample_number(ctid,tcover,header),t_mean(ctid,tcover,header),t_sigma(ctid,tcover,header),t_skewness(ctid,tcover,header));
            ctid=aln->core.tid;
            tcover =(int*)realloc(tcover,sizeof(int)*(*(header->target_len+ctid)));//realloc memory
            for(int i=0;i<(*(header->target_len+ctid));i++)//initialization
                *(tcover+i)=0;
        }
        if(aln->core.l_qseq==28||aln->core.l_qseq==29){
            *(tcover+aln->core.pos+12)+=1;
        }
    }
    if(flag==2) printf("Finished successfully!\n");
    if(flag==0) printf("Read error, please check your bam file.\n");
    bam_destroy1(aln);
    bam_hdr_destroy(header);
    //close file
    sam_close(in);
    fclose(wf);
    //fclose(wt);
    fclose(wc);
    return 0;
}
void write_wig(int tid,int* cover,bam_hdr_t *header,FILE* wf){
    fprintf(wf,"#bedgraph section %s:0-%d\n",*(header->target_name+tid),*(header->target_len+tid));
    int f=0,l=1;//record the first and last position of fragments with the same value
    int n=*(cover);//record the current value
    for(int i=1;i<*(header->target_len+tid)-1;i++){
        if(*(cover+i)==n){
            l++;
        }
        else{
            fprintf(wf,"%s\t%d\t%d\t%d\n",*(header->target_name+tid),f,l,n);
            f=i;
            n=*(cover+i);
            l=i+1;
        }
    }
    if(*(cover+*(header->target_len+tid)-1)==n){
            l++;
    }
    else{
        f=*(header->target_len+tid)-1;
        n=*(cover+*(header->target_len+tid)-1);
        l=*(header->target_len+tid);
    }
    fprintf(wf,"%s\t%d\t%d\t%d\n",*(header->target_name+tid),f,l,n);
}
void write_cover_tsv(int tid,int* cover,bam_hdr_t *header,FILE* wf){//output percentage coverage to tsv file
    fprintf(wf,"%s",*(header->target_name+tid));
    int length=(int)*(header->target_len+tid);
    int bin_size=length/100;
    int middle=length%100;
    for(int i=0;i<middle;i++){
        double sum=0;
        for(int j=0;j<bin_size+1;j++){
            sum+=*(cover+i*(bin_size+1)+j);
        }
        sum=(double)sum/(double)(bin_size+1);
        fprintf(wf,"\t%f",sum);
    }
    for(int i=middle;i<100;i++){
        double sum=0;
        for(int j=0;j<bin_size;j++){
            sum+=*(cover+middle+i*bin_size+j);
        }
        sum=(double)sum/(double)(bin_size);
        fprintf(wf,"\t%f",sum);
    }
    fprintf(wf,"\n");
}
double sample_number(int tid,int* cover,bam_hdr_t*header){
    double cnt=0;
    for(int i=0;i<*(header->target_len+tid);i++)
        cnt+=(*(cover+i));
    return cnt;
}
double t_mean(int tid,int* cover,bam_hdr_t*header){
    double cnt=sample_number(tid,cover,header);
    double mean=0;
    for(int i=0;i<*(header->target_len+tid);i++)
        mean+=(*(cover+i))*(i+1);
    mean=(double)mean/cnt;
    return(mean);
}
double t_sigma(int tid,int* cover,bam_hdr_t*header){
    double sigma=0;
    double cnt=sample_number(tid,cover,header);
    double mean=t_mean(tid,cover,header);
    for(int i=0;i<*(header->target_len+tid);i++)
        sigma+=(*(cover+i))*pow((i+1-mean),2);
    sigma=sigma/(cnt-1);
    sigma=sqrt(sigma);
    return sigma;
}
double t_skewness(int tid,int* cover,bam_hdr_t*header){
    double skewness=0;
    double cnt=sample_number(tid,cover,header);
    double mean=t_mean(tid,cover,header);
    double sigma=t_sigma(tid,cover,header);
    for(int i=0;i<*(header->target_len+tid);i++)
        skewness+=(*(cover+i))*pow((i+1-mean)/sigma,3);
    skewness=skewness/cnt;
    return skewness;
}