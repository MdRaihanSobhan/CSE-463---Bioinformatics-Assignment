#include<bits/stdc++.h>
using namespace std;

string find_random_motif(string str, int k, int n)
{
    int pos = rand() % (n - k);
    return str.substr(pos, k);
}

vector<string> find_profile_most_motifs(int t, int n, int k, string dna[], double profile[][100])
{
    double score, best_score;
    string temp;
    vector<string> profile_most_motifs;

    for(int i = 0; i < t; ++i)
    {
        best_score = 0;
        for(int j = 0; j < n-k; ++j)
        {
            score = 1;
            for(int pos = 0; pos < k; ++pos)
            {
                switch(dna[i][j+pos])
                {
                    case 'A':
                    case 'a': score *= profile[0][pos];
                            break;
                    case 'C':
                    case 'c': score *= profile[1][pos];
                            break;
                    case 'G':
                    case 'g': score *= profile[2][pos];
                            break;
                    case 'T':
                    case 't': score *= profile[3][pos];
                            break;
                }
            }
            if(score > best_score)
            {
                best_score = score;
                temp = dna[i].substr(j, k);
            }
        }
        profile_most_motifs.push_back(temp);
    }

    return profile_most_motifs;
}

void random_motif_search(int n, int t, int k, int it, string dna[])
{
    double profile[5][100];
    vector<string> profile_most_motifs, random_motifs;

    for(int i = 0; i < it; ++i)
    {
        profile_most_motifs.clear();
        random_motifs.clear();

        for(int j = 0; j < t; ++j)
            random_motifs.push_back( find_random_motif(dna[j], k, n) );

        int aa, cc, gg, tt;
        for(int j = 0; j < k; ++j)
        {
            aa = cc = gg = tt = 0;
            for(int l = 0; l < t; ++l)
            {
                switch(random_motifs[l][j])
                {
                    case 'A':
                    case 'a': ++aa;
                            break;
                    case 'C':
                    case 'c': ++cc;
                            break;
                    case 'G':
                    case 'g': ++gg;
                            break;
                    case 'T':
                    case 't': ++tt;
                            break;
                }
            }
            profile[0][j] = aa / (double) t;
            profile[1][j] = cc / (double) t;
            profile[2][j] = gg / (double) t;
            profile[3][j] = tt / (double) t;
        }

        profile_most_motifs = find_profile_most_motifs(t, n, k, dna, profile);
    }

    cout << "Random Motifs:" << endl;
    for(int i = 0; i < t; ++i)
        cout << random_motifs[i] << endl;
    cout << endl;

    char nuces[] = {'A', 'C', 'G', 'T'};
    cout << "Profile:" << endl;
    for(int i = 0; i < 4; ++i)
    {
        cout << nuces[i] << ": ";
        for(int j = 0; j < k; ++j)
            cout << fixed << setprecision(4) << profile[i][j] << " ";
        cout << endl;
    }
    cout << endl;

    cout << "Profile Most Motifs:" << endl;
    for(int i = 0; i < t; ++i)
        cout << profile_most_motifs[i] << endl;

    return;
}

int main()
{
    srand(time(NULL));
    int n, t, k, it;
    string dna[100];
    freopen("RMS_yst08r.txt","r",stdin);
    cin >> n >> t >> k >> it;
    for(int i = 0; i < t; ++i)
        cin >> dna[i];

    random_motif_search(n, t, k, it, dna);

    return 0;
}
/*
10 5 4 1
ttACCTtaac
gATGTctgtc
ccgGCGTtag
cactaACGAg
cgtcagAGGT
*/