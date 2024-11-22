#include <stdio.h>
#include <stdlib.h>
#include <string.h>
struct trans_structure{
    char gene_id[25];
    char transcript_id[25];
    char strand;
    int length;
    int start_coord;//first bp of start codon
    int stop_coord;//last bp of stop codon
    };
struct count_feature{
    char gene_id[25];
    char transcript_id[25];
    char strand;
    int length;
    int UTR5;
    int UTR3;
    int CDS;
    };
void read_wig(FILE *f, int n, int *coverage);
void read_annotation(FILE *f, struct trans_structure *annotation);
int read_head(char*target,int *length,FILE*f);
int length(char*a);
void copystr(char*a,char*b);
void write_coverage(FILE *f, struct count_feature *count,int length);

int main(int argc,char**argv){
    FILE*wig,*tsv,*out1,*out2;//gtf file
    wig=fopen(argv[2],"r");//open wiggle file for scan
    tsv=fopen(argv[1],"r");//open annotation
    out1=fopen(argv[3],"w+");//open output tsv1
    out2=fopen(argv[4],"w+");//open output tsv2
    int length =22500;
    int translen=0;
    char ENST[25];
    struct trans_structure *annotation=(struct trans_structure*)malloc(sizeof(struct trans_structure)*length);
    struct count_feature *count=(struct count_feature*)malloc(sizeof(struct count_feature)*length);
    int a=0;
    read_annotation(tsv,annotation);
    int flag = 0;//record the length of features
    int b=-1;
    fprintf(out2,"transcript_id");
    for(int i=0;i<200;i++){
        fprintf(out2,"\tP%d",i-100);
    }
    fprintf(out2,"\n");
    while(1){
        a=read_head(ENST,&translen,wig);
        if (a==1){
            break;
        }
        if (b==1){
            break;
        }
        int *coverage, outcoverage[200];
        coverage=(int *)malloc(sizeof(int)*translen);
        read_wig(wig,translen,coverage);
        for(int i=0;i<length;i++){
            if(strcmp(annotation[i].transcript_id,ENST)==0){
                copystr(count[flag].gene_id,annotation[i].gene_id);
                copystr(count[flag].transcript_id,annotation[i].transcript_id);
                count[flag].strand=annotation[i].strand;
                count[flag].length=annotation[i].length;
                count[flag].UTR5=0;
                count[flag].UTR3=0;
                count[flag].CDS=0;
                for(int j=0;j<annotation[i].start_coord;j++){
                    count[flag].UTR5+=coverage[j];
                }
                for(int j=annotation[i].stop_coord;j<translen;j++){
                    count[flag].UTR3+=coverage[j];
                }
                for(int j=annotation[i].start_coord;j<annotation[i].stop_coord;j++){
                    count[flag].CDS+=coverage[j];
                }
                flag++;
                for(int j=0;j<200;j++){
                    if((annotation[i].start_coord-101+j)<0)
                        outcoverage[j]=0;
                    else if((annotation[i].start_coord-101+j)>=translen)
                        outcoverage[j]=0;
                    else
                        outcoverage[j]=coverage[annotation[i].start_coord-101+j];
                }
                fprintf(out2,"%s",ENST);
                for(int i=0;i<200;i++){
                    fprintf(out2,"\t%d",outcoverage[i]);
                 }
                fprintf(out2,"\n");
                break;
            }
            
        }
        
        free(coverage);
    }
    write_coverage(out1,count,flag);  
    fclose(wig);
    fclose(tsv);
    fclose(out1);
    fclose(out2);
}
void read_annotation(FILE *f, struct trans_structure *annotation){
    int a=fscanf(f,"%s\t%s\t%c\t%d\t%d\t%d\n",annotation->gene_id,annotation->transcript_id,&(annotation->strand),&(annotation->length),&(annotation->start_coord),&(annotation->stop_coord));
    while(a!=EOF){
        a=fscanf(f,"%s\t%s\t%c\t%d\t%d\t%d\n",annotation->gene_id,annotation->transcript_id,&(annotation->strand),&(annotation->length),&(annotation->start_coord),&(annotation->stop_coord));
        annotation++;
    }
}
void read_wig(FILE *f, int n, int *coverage){
    int start,end,value;
        char ENSEMBL[25];
    while(1){
        fscanf(f,"%s\t%d\t%d\t%d\n",ENSEMBL,&start,&end,&value);
        for(int i=start;i<end;i++)
            *(coverage+i)=value;
        if(end==n)
            break;
    }
}
int read_head(char*target,int *length,FILE*f){
    int cnt=0;
    for(int i=0;i<18;i++){
        int a=fgetc(f);
        if(a==EOF)
            return 1;
    }
    while(1){
        char c=fgetc(f);
        if(c==':'){
            *(target+cnt)='\0';
            break;
        }
        *(target+cnt)=c;
        cnt++;
    }
    for(int i=0;i<2;i++)
        fgetc(f);
    fscanf(f,"%d\n",length);
    return 0;
}
int length(char*a){
    int cnt=0;
    while(a[cnt]!='\0'){
        if(a[cnt]=='\n'){
            a[cnt]='\0';
            break;
        }
        cnt++;
    }
    return cnt;
}
void copystr(char*a,char*b){//b--->a
    int la=length(a),lb=length(b);
    if(la>lb){
        for(int i=0;i<la;i++){
            if(i<lb)
                *(a+i)=*(b+i);
            else
                *(a+i)='\0';
        }
    }
    else{
        for(int i=0;i<lb;i++){
            *(a+i)=*(b+i);
        *(a+lb)='\0';
    }
    }
}
void write_coverage(FILE *f, struct count_feature *count,int length){
    fprintf(f,"gene_id\ttranscript_id\tstrand\tlength\tUTR5\tCDS\tUTR3\n");
    for(int i=0;i<length;i++){
        fprintf(f,"%s\t%s\t%c\t%d\t%d\t%d\t%d\n",count[i].gene_id,count[i].transcript_id,count[i].strand,count[i].length,count[i].UTR5,count[i].CDS,count[i].UTR3);
    }
}
