#pragma warning(disable:4996)
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <iostream>
#include<complex>

# define complex std::complex<double>
# define cycles 1
# define DELTA_T 0

using namespace std;

#define S 262144                   //Cache size
#define MX_LS 65536            //MAX_LINES
#define MX_WS 16                //MAX_WAYS
#define data_bus_width 32   // 32 bit ------ 4 byte

int LRU[MX_LS][MX_WS];       
int TAG[MX_LS][MX_WS];       //TAG
bool v[MX_LS][MX_WS];        //Valid bit
bool d[MX_LS][MX_WS];        //Dirty bit

uint8_t N;                          //No of ways
int BL;                             //Burst Length
uint8_t block_size_B;               //Block size
uint32_t L, block_bits;             //no of lines , block bits

uint32_t tag, tag_bits;             //tag
uint32_t line, Line, line_bits;    //line

uint32_t RBC = 0;
uint32_t RMC = 0;
uint32_t RBRC = 0;
uint32_t RBRDC = 0;
uint32_t RBHC = 0;
uint32_t RBMC = 0;         //read counters
uint32_t WTC = 0;
uint32_t WBC = 0;
uint32_t WMC = 0;
uint32_t WBHC = 0;
uint32_t WBMC = 0;
uint32_t WBRC = 0;
uint32_t WBRDC = 0;        //write counters
uint32_t FBC = 0;;         //flush counters



FILE* to_print;

bool H_M, wtb, wtna, wta;       

uint8_t WB = 0;
uint8_t WTNA = 0;
uint8_t WTA = 0;

void parameter_conversion(int BL, uint8_t N);

bool cache_miss_or_hit(uint32_t tag, uint32_t line);

uint32_t calculate_tag_bits(uint32_t address);

uint32_t calculate_line_bits(uint32_t address);

uint8_t LRU_index(uint32_t line);
void to_update_LRU(uint32_t i, uint32_t line);
void Read_Block(uint32_t tag, uint32_t line);
void Read_Memory(uint32_t* address, uint32_t bytes);
void Write_Block(uint32_t tag, uint32_t line);
void write_memory(uint32_t* address, uint32_t bytes);
void clear_cache();
void clear_counters();
void Cache_Flush();
void print_counters();
void clear_string(char* string);
void csv_file_out();
void Radix2FFT(complex data[], int nPoints, int nBits);




int main() 
{
    uint32_t i;
    uint32_t points = 32768; 
    complex data[32768];

    uint32_t x[32768];
    printf("FFT TEST APP\r\n");

     to_print = fopen("Counters.csv", "w+");
    fprintf(to_print, "BL,N,WS,RMC,RBC,RBHC,RBMC,RBRC,RBRDC,WMC,WBC,WBHC,WBMC,WTC,WBRC,WBRDC,FBC \n");
   

    wtb = false;
     wtna = false;
    wta = false;

    for (BL = 1; BL <= 8; BL += BL)
    {
        for (N = 1; N <= 16; N += N)
        {
           for (int write_strtaegy = 1; write_strtaegy <= 3; write_strtaegy++)
            {
                parameter_conversion(BL, N);
                if (write_strtaegy == 1)
                {
                    wtb = true;
                }
                else if (write_strtaegy == 2)
                {
                    wtb = false;
                    wta = true;
                }
                else
                {
                    wtb = false;
                    wta = false;
                    wtna = true;
                }
               
                clear_cache();
                clear_counters();
               
                               
                 write_memory(&i, sizeof(i));
                 for (i = 0; i < points; i++)
                 {
                     data[i] = complex(cos(2.0 * 3.1416 * float(i) / float(points) * cycles), 0.0);
                     x[i] = data[i].real(); // taking the real values from data only.

                     Read_Memory(&i, sizeof(i));
                     write_memory(&x[i], sizeof(i));
                     write_memory(&i, sizeof(i));

                     for (i = 0; i < points; i++)
                     {
                         data[i] = complex(0.0, 0.0);
                         data[DELTA_T] = complex(1.0, 0.0);
                         Read_Memory(&i, sizeof(i));
                         write_memory(&x[i], sizeof(i));
                         write_memory(&i, sizeof(i)); 

                     }
                         int bits = ceil(log(points) / log(2));
                         Radix2FFT(data, points, bits);
                         write_memory(&i, sizeof(i));
                         for (i = 0; i < points; i++)
                         {
                             Read_Memory(&i, sizeof(i));
                             write_memory(&x[i], sizeof(i));
                             write_memory(&i, sizeof(i));
                             if (data[i].imag() >= 0.0)
                             {
                                 printf("x[%d] = %2.4lf + j%2.4lf\n", i, data[i].real(), data[i].imag());
                             }
                             else
                             {
                                 printf("x[%d] = %2.4lf - j%2.4lf\n", i, data[i].real(), -data[i].imag());
                             }
                         }                 
                  
                 }
                 Cache_Flush();
                 print_counters();
                 csv_file_out();
            }

        }

    }

    fclose(to_print);
    return 0;
}




void parameter_conversion(int BL, uint8_t N)
{
    block_size_B = BL * (data_bus_width / 8);
    block_bits = (uint8_t)log2(block_size_B);

    L = S / (4 * N * BL);
    line_bits = (uint32_t)log2(L);

    tag_bits = data_bus_width - block_bits - line_bits;
}

//Cache hit or miss logic
bool cache_miss_or_hit(uint32_t tag, uint32_t line) {

    H_M = false;
    uint32_t i, j;

    for (i = 0; i < N; i++) {
        j = TAG[i][line];
        if (j == tag) {          //if tag matches the given tag at that particular line
            H_M = true;
            break;
        }
        else
        {
            H_M = false;
        }
    }
    return H_M;
}

//calculate tag from  address
uint32_t calculate_tag_bits(uint32_t address)
{
    uint32_t cal_tag;
    uint8_t x = block_bits + line_bits;
    cal_tag = address >> x;
    return cal_tag;
}


//calculate line from address
uint32_t calculate_line_bits(uint32_t address)
{

    uint32_t cal_line;
    uint8_t x = block_bits + tag_bits;
    cal_line = (address << tag_bits) >> x;
    return cal_line;
}

//LRU Logic 
uint8_t LRU_index(uint32_t line) 
{

    uint8_t i, z, x = 0;
    for (i = 0; i < N; i++) 
    {
        z = LRU[i][line];
        if (z == (N - 1))          //finding the LRU value at particular line
        {
            x = i;
            break;
        }
    }
    return x;
}

//Update LRU Logic
void to_update_LRU(uint32_t i, uint32_t line) 
{
    uint32_t x = LRU[i][line];
    LRU[i][line] = 0;
    for (int m = 0; m < N; m++)

    {
        if (LRU[m][line] < x)
        {
            LRU[m][line] = LRU[m][line] + 1;
        }
    }

}

void Read_Block(uint32_t tag, uint32_t line)
{
    uint8_t LRU_index_rd = 0;
    RBC++;
    H_M = cache_miss_or_hit(tag, line);
    LRU_index_rd = LRU_index(line);
    if (H_M == true)
    {
        RBHC++;
    }
    else
    {
        RBMC++;
        LRU_index_rd = LRU_index(line);
        if (v[LRU_index_rd][line] == 1)
        {
            RBRC++;

            if (d[LRU_index_rd][line] == 1)
            {
                RBRDC++;
                v[LRU_index_rd][line] = 0;
                d[LRU_index_rd][line] = 0;
            }
        }
        TAG[LRU_index_rd][line] = tag;
        v[LRU_index_rd][line] = 1;
        d[LRU_index_rd][line] = 0;
    }
    to_update_LRU(LRU_index_rd, line);

}

//Read Memory logic
void Read_Memory(uint32_t* address, uint32_t bytes)
{

    RMC++;
    uint32_t New_address = (uint32_t)address;
    int old_line = -1;
    for (int i = 0; i < bytes; i++)
    {
        tag = calculate_tag_bits(New_address);
        line = calculate_line_bits(New_address);
        if (line != old_line) {
            old_line = line;
            Read_Block(tag, line);
        }
        New_address++;
    }

}

//write block 
void Write_Block(uint32_t tag, uint32_t line)
{

    WBC++;
    uint8_t LRU_index_wt = 0;
    H_M = cache_miss_or_hit(tag, line);
    if (H_M == false)
    {
        WBMC++;
        LRU_index_wt = LRU_index(line);
        if (v[LRU_index_wt][line] == 0 && ((wtb) || (wta))) // if not free 
        {
            WBRC++;

            if (d[LRU_index_wt][line] == 1)
            {
                WBRDC++;
                d[LRU_index_wt][line] = 0;
            }

            TAG[LRU_index_wt][line] = tag; // tag update
            v[LRU_index_wt][line] = 1;
            to_update_LRU(LRU_index_wt, line); // update lru

        }
    }
    else {
        WBHC++;
        v[LRU_index_wt][line] = 1;
        to_update_LRU(LRU_index_wt, line);
    }
    if (wtb)
    {
        v[LRU_index_wt][line] = 1;
        d[LRU_index_wt][line] = 1;
    }
}


void write_memory(uint32_t* address, uint32_t bytes)
{

    WMC++;
    uint32_t New_address = (uint32_t)address;
    int old_line = -1;

    for (int i = 0; i < bytes; i++)
    {
        tag = calculate_tag_bits(New_address);
        line = calculate_line_bits(New_address);
        if (line != old_line) {
            old_line = line;
            Write_Block(tag, line);
        }
        New_address++;
    }
    if ((wta) || (wtna))
    {

        WTC++;
    }
}

// cache clear logic
void clear_cache()
{
    int i, j;

    for (i = 0; i < MX_WS; i++)
    {
        for (j = 0; j < MX_LS; j++)
        {
            TAG[i][j] = 0;
            LRU[i][j] = 0;
            v[i][j] = 0;
            d[i][j] = 0;
        }
    }
}


//clear  all thecounters
void clear_counters()
{

    RMC = 0;
    RBC = 0;
    RBRC = 0;
    RBRDC = 0;
    RBHC = 0;
    RBMC = 0;
    WMC = 0;
    WTC = 0;
    WBC = 0;
    WBHC = 0;
    WBMC = 0;
    WBRC = 0;
    WBRDC = 0;
    FBC = 0;

}

//cache flush logic
void Cache_Flush()
{
    bool ready = false;
    int i, j;
    if (!wtna)
    {
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < L; j++)
            {
                if (v[i][j] == 1 || d[i][j] == 1)
                {
                    ready = true;
                }
                FBC++;
            }

        }
    }
}

void print_counters()
{

    printf("  \n");
    printf("BL %d \n", BL);
    printf("N %d \n", N);
    if (wtb) 
    {
        printf("Write-back \n");
    }
    else if (wtna) 
    {
        printf("Write Through Non Allocate \n");
    }
    else {
        printf("Write Through Allocate \n");
    }
    printf("RMC %d \n", RMC);
    printf("RBC %d \n", RBC);
    printf("RBHC %d \n", RBHC);
    printf("RBMC %d \n", RBMC);
    printf("RBRDC %d \n", RBRC);
    printf("RBDC %d \n", RBRDC);
    printf("WMC %d \n", WMC);
    printf("WTC %d \n", WTC);
    printf("WBC %d \n", WBC);
    printf("WBHC %d \n", WBHC);
    printf("WBMC %d \n", WBMC);
    printf("WBRC %d \n", WBRC);
    printf("WBRDC %d \n", WBRDC);
    printf("FBC %d \n", FBC);
    printf("  \n");
}
void clear_string(char* string)
{
    int i = 0;
    for (i = 0; i < strlen(string); i++)
    {
        string[i] = '\0';
    }
}

void csv_file_out()
{

    fprintf(to_print, "%d,", BL);
    fprintf(to_print, "%d,", N);

    if (wtb) {
        char string1[3] = ("WB");
        fprintf(to_print, "%s,", string1);
        clear_string(string1);
    }
    else if (wta) {
        char string2[4] = ("WTA");
        fprintf(to_print, "%s,", string2);
        clear_string(string2);
    }
    else {
        char string3[5] = ("WTNA");
        fprintf(to_print, "%s,", string3);
        clear_string(string3);
    }

    fprintf(to_print, "%d,", RMC);
    fprintf(to_print, "%d,", RBC);
    fprintf(to_print, "%d,", RBHC);
    fprintf(to_print, "%d,", RBMC);
    fprintf(to_print, "%d,", RBRC);
    fprintf(to_print, "%d,", RBRDC);
    fprintf(to_print, "%d,", WBC);
    fprintf(to_print, "%d,", WMC);
    fprintf(to_print, "%d,", WBHC);
    fprintf(to_print, "%d,", WBMC);
    fprintf(to_print, "%d,", WTC);
    fprintf(to_print, "%d,", WBRC);
    fprintf(to_print, "%d,", WBRDC);
    fprintf(to_print, "%d,\n", FBC);
}


void Radix2FFT(complex data[], int nPoints, int nBits)
{
    // cooley-tukey radix-2, twiddle factor
    // adapted from Fortran code from Burrus, 1983
#pragma warning (disable: 4270)
    uint32_t i, j, k, l;
    int nPoints1, nPoints2;
    complex cTemp, cTemp2;
    nPoints2 = nPoints;
    uint32_t y[32768];
    uint32_t z[32768];
    write_memory(&k, sizeof(k));
    for (k = 1; k <= nBits; k++)
    {
        y[k] = data[k].real();
        Read_Memory(&k, sizeof(k));
        write_memory(&y[k], sizeof(k));
        write_memory(&k, sizeof(k));

        nPoints1 = nPoints2;
        nPoints2 /= 2;
        // Compute differential angles
        double dTheta = 2 * 3.14159257 / nPoints1;
        double dDeltaCos = cos(dTheta);
        double dDeltaSin = sin(dTheta);
        // Initialize angles
        double dCos = 1;
        double dSin = 0;
        // Perform in-place FFT
        write_memory(&j, sizeof(j));
        for (j = 0; j < nPoints2; j++)
        {
            Read_Memory(&j, sizeof(j));
            write_memory(&y[j], sizeof(j));
            write_memory(&j, sizeof(j));
            i = j;
            write_memory(&i, sizeof(i));
            while (i < nPoints)
            {
                Read_Memory(&i, sizeof(i));
                write_memory(&y[i], sizeof(i));
                write_memory(&i, sizeof(i));

                l = i + nPoints2;
                cTemp = data[i] - data[l];
                cTemp2 = data[i] + data[l];
                data[i] = cTemp2;
                z[i] = data[i].real();
                Read_Memory(&i, sizeof(i));
                write_memory(&y[i], sizeof(i));
                write_memory(&i, sizeof(i));

                cTemp2 = complex(dCos * cTemp.real() + dSin * cTemp.imag(), dCos * cTemp.imag() - dSin * cTemp.real());
                data[l] = cTemp2;
                Read_Memory(&l, sizeof(l));
                write_memory(&y[l], sizeof(l));
                write_memory(&l, sizeof(l));

                i += nPoints1;
            }
            double dTemp = dCos;
            dCos = dCos * dDeltaCos - dSin * dDeltaSin;
            dSin = dTemp * dDeltaSin + dSin * dDeltaCos;
        }
    }
    // Convert Bit Reverse Order to Normal Ordering
    j = 0;
    nPoints1 = nPoints - 1;
    write_memory(&i, sizeof(i));
    for (i = 0; i < nPoints1; i++)
    {
        Read_Memory(&i, sizeof(i));
        write_memory(&y[i], sizeof(i));
        write_memory(&i, sizeof(i));
        if (i < j)
        {
            Read_Memory(&j, sizeof(j));
            write_memory(&y[j], sizeof(j));
            write_memory(&j, sizeof(j));

            cTemp = data[j];

            Read_Memory(&j, sizeof(j));
            write_memory(&y[j], sizeof(j)); // check once
            write_memory(&j, sizeof(j));

            cTemp2 = data[i];
            Read_Memory(&i, sizeof(i));
            write_memory(&y[i], sizeof(i)); // check once
            write_memory(&i, sizeof(i));
            data[i] = cTemp;
            data[j] = cTemp2;
        }
        k = nPoints / 2;
        write_memory(&k, sizeof(k));
        while (k <= j)
        {
            Read_Memory(&k, sizeof(k));
            write_memory(&y[k], sizeof(k)); 
            write_memory(&k, sizeof(k));

            j -= k;
            k /= 2;
        }
        j += k;
    }
#pragma warning(default: 4270)
}
