#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <fstream>

using namespace std;

class Algorithm {
    public:
    string name;
    int rank;
    double fastest = 100000;
    double slowest = 0;
    double current = 0;
    int id;
    int matches;
    int points;
    int wins;
    int losses;
    int ties;
    double totalTime;

    Algorithm(string n, int i) {
        name = n;
        id = i;
    }

    virtual void sort(vector<int> a){
        return;
    }

    void display(){
        cout << "Name: " << name << "\nFastest Time: " << fastest / 1000 << "s" << "\nSlowest Time: " << slowest / 1000 << "s" << "\nCurrent Time: " << current / 1000 << "s" << endl;
    }
};

class BubbleSort : public Algorithm {
    public:
    BubbleSort() : Algorithm("Bubble Sort", 2) {};

    void sort(vector<int> a){
        auto start = chrono::high_resolution_clock::now();

        //Algorithm: 
        int size = a.size();
        for (int i = 0; i < size - 1; i++){
            for(int j = 0; j < size - i - 1; j++){
                if(a[j] > a[j + 1]){
                    int temp = a[j];
                    a[j] = a[j + 1];
                    a[j + 1] = temp;
                }
            }
        }


        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration = end - start;
        double time = duration.count();
        current = time;
        if (time < fastest) fastest = time;
        if (time > slowest) slowest = time;
    }
};

class SelectionSort : public Algorithm {
    public:
    SelectionSort() : Algorithm("Selection Sort", 1) {};

    void sort(vector<int> a){
        auto start = chrono::high_resolution_clock::now();

        //Algorithm: 
        int size = a.size();
        for(int i = 0; i < size - 1; i++){
            int minIndex = i;
            for(int j = i + 1; j < size; j++){
                if(a[j] < a[minIndex]){
                    minIndex = j;
                }
            }
            if(minIndex != i){
                int temp = a[i];
                a[i] = a[minIndex];
                a[minIndex] = temp;
            }
        }

        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration = end - start;
        double time = duration.count();
        current = time;
        if (time < fastest) fastest = time;
        if (time > slowest) slowest = time;
    }
};

class MergeSort : public Algorithm {
    private:
    void merge(vector<int>& arr, int left, int mid, int right){
        int n1 = mid - left + 1;
        int n2 = right - mid;
        vector<int> leftArr(n1);
        vector<int> rightArr(n2);

        for(int i = 0; i < n1; i++) leftArr[i] = arr[left + i];
        for(int j = 0; j < n2; j++) rightArr[j] = arr[mid + 1 + j];

        int i = 0, j = 0, k = left;
        while(i < n1 && j < n2) {
            if(leftArr[i] <= rightArr[j]){
                arr[k] = leftArr[i];
                i++;
            } else {
                arr[j] = rightArr[j];
                j++;
            }
            k++;
        }

        while (i < n1){
            arr[k] = leftArr[i];
            i++;
            k++;
        }

        while (j < n2){
            arr[k] = rightArr[j];
            j++;
            k++;
        }
    }

    void merge_sort(vector<int>& arr, int left, int right) {
        if(left < right){
            int mid = left + (right - left) / 2;

            merge_sort(arr, left, mid);
            merge_sort(arr, mid + 1, right);
            merge(arr, left, mid, right);
        }
    }

    public:
    MergeSort() : Algorithm("Merge Sort", 3) {};

    void sort(vector<int> a){
        auto start = chrono::high_resolution_clock::now();

        //Algorithm: 
        merge_sort(a, 0, a.size() - 1);

        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration = end - start;
        double time = duration.count();
        current = time;
        if (time < fastest) fastest = time;
        if (time > slowest) slowest = time;
    }
};

class InsertionSort : public Algorithm {
    public:
    InsertionSort() : Algorithm("Insertion Sort", 4) {};

    void sort(vector<int> a){
        auto start = chrono::high_resolution_clock::now();

        //Algorithm: 
        int size = a.size();
        for (int i = 1; i < size; i++){
            int key = a[i];
            int j = i - 1;

            while(j >= 0 && a[j] > key){
                a[j + 1] = a[j];
                j--;
            }
            a[j + 1] = key;
        }

        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration = end - start;
        double time = duration.count();
        current = time;
        if (time < fastest) fastest = time;
        if (time > slowest) slowest = time;
    }
};
class CycleSort : public Algorithm {
    public:
    CycleSort() : Algorithm("Cycle Sort", 5) {};

    void sort(vector<int> a){
        auto start = chrono::high_resolution_clock::now();
        int n = a.size();
        int writes = 0;
    
        for (int cycle_start = 0; cycle_start <= n - 2; cycle_start++) {
            int item = a[cycle_start];
            int pos = cycle_start;
            for (int i = cycle_start + 1; i < n; i++)
                if (a[i] < item)
                    pos++;
            if (pos == cycle_start)
                continue;
            while (item == a[pos])
                pos += 1;

            if (pos != cycle_start) {
                swap(item, a[pos]);
                writes++;
            }

            while (pos != cycle_start) {
                pos = cycle_start;
    
                for (int i = cycle_start + 1; i < n; i++)
                    if (a[i] < item)
                        pos += 1;
    
                while (item == a[pos])
                    pos += 1;
    
                if (item != a[pos]) {
                    swap(item, a[pos]);
                    writes++;
                }
            }
        }

        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration = end - start;
        double time = duration.count();
        current = time;
        if (time < fastest) fastest = time;
        if (time > slowest) slowest = time;
    }
};

class QuickSort : public Algorithm {
    private:
    int partition(vector<int>& arr, int low, int high) {
        int pivot = arr[high];
        int i = low - 1;
        for (int j = low; j <= high - 1; j++) {
            if (arr[j] < pivot) {
                i++;
                swap(arr[i], arr[j]);
            }
        }
        swap(arr[i + 1], arr[high]);  
        return i + 1;
    }

    void quickSort(vector<int>& arr, int low, int high) { 
        if (low < high) {
            int pi = partition(arr, low, high);
            quickSort(arr, low, pi - 1);
            quickSort(arr, pi + 1, high);
        }
    }

    public:
    QuickSort() : Algorithm("Quick Sort", 6) {};


    void sort(vector<int> a){
        auto start = chrono::high_resolution_clock::now();

        //Algorithm: 
        int size = a.size();
        quickSort(a, 0, size - 1);

        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration = end - start;
        double time = duration.count();
        current = time;
        if (time < fastest) fastest = time;
        if (time > slowest) slowest = time;
    }
};

class HeapSort : public Algorithm {
    private:
    void heapify(vector<int>& arr, int n, int i){
        int largest = i;
        int l = 2 * i + 1;
        int r = 2 * i + 2;
        if (l < n && arr[l] > arr[largest])
            largest = l;
        if (r < n && arr[r] > arr[largest])
            largest = r;
        if (largest != i) {
            swap(arr[i], arr[largest]);
            heapify(arr, n, largest);
        }
    }

    void heapSort(vector<int>& arr){
        int n = arr.size();
        for (int i = n / 2 - 1; i >= 0; i--)
            heapify(arr, n, i);
        for (int i = n - 1; i > 0; i--) {
            swap(arr[0], arr[i]);
            heapify(arr, i, 0);
        }
    }


    public:
    HeapSort() : Algorithm("Heap Sort", 7) {};


    void sort(vector<int> a){
        auto start = chrono::high_resolution_clock::now();

        //Algorithm: 
        int size = a.size();
        heapSort(a);

        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration = end - start;
        double time = duration.count();
        current = time;
        if (time < fastest) fastest = time;
        if (time > slowest) slowest = time;
    }
};
class ThreeWayMergeSort : public Algorithm {
    private:
    void merge(vector<int> &gArray, int low, int mid1, int mid2, int high, vector<int> &destArray) 
    { 
        int i = low, j = mid1, k = mid2, l = low; 
        while ((i < mid1) && (j < mid2) && (k < high)) { 
            if(gArray[i] < gArray[j]){
                if(gArray[i] < gArray[k]){
                    destArray[l++] = gArray[i++];
                }
                else{
                    destArray[l++] = gArray[k++];
                }
            }
            else{
                if(gArray[j] < gArray[k]){
                    destArray[l++] = gArray[j++];
                }
                else{
                    destArray[l++] = gArray[k++];
                }
            }
        } 

        while ((i < mid1) && (j < mid2)) { 
            if(gArray[i] < gArray[j]){
                destArray[l++] = gArray[i++];
            }
            else {
                destArray[l++] = gArray[j++];
            }
        } 
        while ((j < mid2) && (k < high)) { 
            if(gArray[j] < gArray[k]){
                destArray[l++] = gArray[j++];
            }
            else{
                destArray[l++] = gArray[k++];
            } 
        } 

        while ((i < mid1) && (k < high)) { 
            if(gArray[i] < gArray[k]){
                destArray[l++] = gArray[i++];
            }
            else{
                destArray[l++] = gArray[k++];
            } 
        } 
        while (i < mid1) 
            destArray[l++] = gArray[i++]; 
    
        while (j < mid2) 
            destArray[l++] = gArray[j++]; 
    
        while (k < high) 
            destArray[l++] = gArray[k++]; 
    } 

    void mergeSort3WayRec(vector<int> &gArray, int low, int high, vector<int> &destArray) { 
        if (high - low < 2) 
            return; 

        int mid1 = low + ((high - low) / 3); 
        int mid2 = low + 2 * ((high - low) / 3) + 1; 

        mergeSort3WayRec(destArray, low, mid1, gArray); 
        mergeSort3WayRec(destArray, mid1, mid2, gArray); 
        mergeSort3WayRec(destArray, mid2, high, gArray); 

        merge(destArray, low, mid1, mid2, high, gArray); 
    }
    
    void mergeSort3Way(vector<int> &gArray, int n) 
    { 

        if (n == 0) 
            return; 
    
        vector<int> fArray(n); 

        for (int i = 0; i < n; i++) 
            fArray[i] = gArray[i]; 

        mergeSort3WayRec(fArray, 0, n, gArray); 

        for (int i = 0; i < n; i++) 
            gArray[i] = fArray[i]; 
    } 


    public:
    ThreeWayMergeSort() : Algorithm("3-Way Merge", 8) {};


    void sort(vector<int> a){
        auto start = chrono::high_resolution_clock::now();

        //Algorithm: 
        int size = a.size();
        mergeSort3Way(a, size);

        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration = end - start;
        double time = duration.count();
        current = time;
        if (time < fastest) fastest = time;
        if (time > slowest) slowest = time;
    }
};
class RadixSort : public Algorithm {
    private:
    int getMax(const vector<int>& arr) {
        if (arr.empty()) return 0;
        int mx = arr[0];
        for (size_t i = 1; i < arr.size(); i++) {
            mx = max(mx, arr[i]);
        }
        return mx;
    }

    void countSort(vector<int>& arr, long exp) {
        size_t n = arr.size();
        if (n == 0) return;

        try {
            vector<int> output(n);
            vector<int> count(10, 0);

            // Store count of occurrences
            for (size_t i = 0; i < n; i++) {
                int digit = (arr[i] / exp) % 10;
                if (digit >= 0 && digit < 10) {
                    count[digit]++;
                }
            }

            // Change count[i] so that count[i] contains actual position
            for (int i = 1; i < 10; i++) {
                count[i] += count[i - 1];
            }

            // Build the output array
            for (long i = n - 1; i >= 0; i--) {
                int digit = (arr[i] / exp) % 10;
                if (digit >= 0 && digit < 10 && count[digit] > 0) {
                    size_t pos = count[digit] - 1;
                    if (pos < n) {
                        output[pos] = arr[i];
                        count[digit]--;
                    }
                }
            }

            // Copy back to original array
            arr = output;
        } catch (const std::exception& e) {
            // Handle any potential exceptions
            return;
        }
    }

    void radixsort(vector<int>& arr) {
        if (arr.empty()) return;

        // Find the maximum number to know number of digits
        int m = getMax(arr);
        if (m <= 0) return;

        // Do counting sort for every digit
        for (long exp = 1; m / exp > 0 && exp <= 1000000000; exp *= 10) {
            countSort(arr, exp);
        }
    }

    public:
    RadixSort() : Algorithm("Radix Sort", 9) {};

    void sort(vector<int> a) {
        auto start = chrono::high_resolution_clock::now();
        
        try {
            // Create a working copy
            vector<int> temp = a;
            radixsort(temp);
        } catch (const std::exception& e) {
            // Handle any potential exceptions
        }
        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration = end - start;
        double time = duration.count();
        current = time;
        if (time < fastest) fastest = time;
        if (time > slowest) slowest = time;
    }
};

class CountingSort : public Algorithm {
    public:
    CountingSort() : Algorithm("Counting Sort", 10) {};

    void sort(vector<int> a) {
        auto start = chrono::high_resolution_clock::now();

        vector<int> temp = a;  // Create working copy
        int N = temp.size();
        
        // Find maximum element
        int M = 0;
        for (int i = 0; i < N; i++)
            M = max(M, temp[i]);
            
        vector<int> countArray(M + 1, 0);
        vector<int> outputArray(N);
        
        // Store count of each element
        for (int i = 0; i < N; i++)
            countArray[temp[i]]++;
            
        // Modify countArray to store actual positions
        for (int i = 1; i <= M; i++)
            countArray[i] += countArray[i - 1];
            
        // Build output array
        for (int i = N - 1; i >= 0; i--) {
            outputArray[countArray[temp[i]] - 1] = temp[i];
            countArray[temp[i]]--;
        }
        
        // Copy back to working array
        temp = outputArray;
        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration = end - start;
        double time = duration.count();
        current = time;
        if (time < fastest) fastest = time;
        if (time > slowest) slowest = time;
    }
};
class Tournament {
    int matchNo = 1;
     vector<Algorithm*> algorithms;

public:
    Tournament() {
        algorithms.push_back(new SelectionSort());
        algorithms.push_back(new RadixSort());
        algorithms.push_back(new BubbleSort());
        algorithms.push_back(new MergeSort());
        algorithms.push_back(new InsertionSort());
        algorithms.push_back(new QuickSort());
        algorithms.push_back(new HeapSort());
        algorithms.push_back(new CycleSort());
        algorithms.push_back(new ThreeWayMergeSort());
        algorithms.push_back(new CountingSort());
    }

    ~Tournament() {
        for (auto algo : algorithms) {
            delete algo; 
        }
    }

    void simluateMatch(Algorithm* team1, Algorithm* team2, vector<int> a){
        team1->matches++;
        team2->matches++;
        team1->sort(a);
        team2->sort(a);

        cout << team1->name << "\tvs\t" << team2->name << endl;
        cout << std::setprecision(7) << team1->current / 1000 << "s" << "\t\t" << std::setprecision(7) << team2->current / 1000 << "s" << endl;

        cout << "WINNER: ";

        if(team1->current < team2->current) {
            cout << team1->name << endl;
            team1->points += 3;
            team1->wins++;
            team2->losses++;
        }
        else if(team2->current < team1->current){
            cout << team2->name << endl;
            team2->points += 3;
            team2->wins++;
            team1->losses++;
        } 
        else if (team1->current == team2->current){
            cout << "TIED" << endl;
            team1->points += 1;
            team2->points += 1;
            team1->ties++;
            team2->ties++;
        }
    }

    void roundRobin(vector<int> a){
        for (size_t i = 0; i < algorithms.size(); i++){
            for(size_t j = 0; j < algorithms.size(); j++){
                if( i != j){
                cout << "MATCH " << matchNo << ": " <<  algorithms[i]->name << " vs " << algorithms[j]->name << endl;
                simluateMatch(algorithms[i], algorithms[j], a);
                cout << endl;
                matchNo++;
                }
            }
            cout << "End of Round " << i + 1 << endl;
            printLeaderboard();
            cout << endl;
        }
    }



    void printLeaderboard() {
        // Create a copy of the algorithms vector
        vector<Algorithm*> sortedAlgorithms = algorithms;

        // Sort the copied vector in descending order based on points
        sort(sortedAlgorithms.begin(), sortedAlgorithms.end(), [](Algorithm* a, Algorithm* b) {
            return a->points > b->points;
        });

        // Display the leaderboard with positions
        std::cout << "\nPoints Table:\n";
    std::cout << "Pos\tAlgo\t\tMatches\tWins\tLosses\tTies\tFastest\t\tSlowest\t\tLast\t\tPoints\n";
    std::cout << "---------------------------------------------------------\n";
    
    int position = 1;
    for (const auto &alg : sortedAlgorithms) {
        std::cout << std::setw(3) << position++ << "\t" 
                  << std::setw(15) << alg->name << "\t" 
                  << std::setw(7) << alg->matches << "\t" 
                  << std::setw(4) << alg->wins << "\t" 
                  << std::setw(6) << alg->losses << "\t" 
                  << std::setw(5) << alg->ties << "\t" 
                  << std::setw(10) << std::fixed << std::setprecision(7) << (alg->fastest / 1000) << "s\t" 
                  << std::setw(10) << std::fixed << std::setprecision(7) << (alg->slowest / 1000) << "s\t" 
                  << std::setw(10) << std::fixed << std::setprecision(7) << (alg->current / 1000) << "s\t" 
                  << alg->points << "\n";
    }
    }

};



void readFile(const string& filename, vector<int>& arr){
    ifstream file(filename);
    if(!file.is_open()){
        return;
    }

    int num;
    while (file >> num){
        arr.push_back(num);
    }
    file.close();
}

int main(void){
    vector<int> arr;



   Tournament tourney;
   arr.clear();
   readFile("sorted.txt", arr);

   tourney.roundRobin(arr);

   cout << "END OF SORTED ROUND" << endl;

   arr.clear();
   readFile("reverse.txt", arr);

   tourney.roundRobin(arr);

   cout << "END OF REVERSE ROUND" << endl;

   arr.clear();
   readFile("input.txt", arr);

   tourney.roundRobin(arr);
   cout << "END OF REVERSE ROUND" << endl;

   tourney.printLeaderboard();
}