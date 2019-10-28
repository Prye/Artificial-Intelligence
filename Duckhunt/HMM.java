import java.io.InputStream;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

public class HMM {
    private double[][] A;   // Transition possibility matrix
    private double[][] B;   // Observation probability matrix
    private double[] pi;    // Initial probability pi
    private Integer N;   // Number of hidden states
    private Integer M;     // Number of observations;
    private Integer[] obs;      // Observations
    private double[][] alpha;   // The param for EM algorithm
    private double[][] beta;    // The param for EM algorithm
    private Stack<Double> scales;
    private int T;
    private int maxIter = 25;

    public double[][] getA() {
        return A;
    }

    public double[][] getB() {
        return B;
    }

    public double[] getPi() {
        return pi;
    }

    public double[][] getAlpha() {
        return alpha;
    }

    public int getT() {
        return T;
    }

    public HMM (int N, int M) {
        this.N = N;
        this.M = M;
        A = new double[N][N];
        B = new double[N][M];
        pi = new double[N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = 1 / N + Math.random() * 0.01 / N;
            }
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                B[i][j] = 1 / N + Math.random() * 0.5 / N;
            }
        }
        for (int i = 0; i < N; i++) {
            pi[i] = 1 / N + Math.random() * 0.01 / N;
        }
        Matrices.normalize(A);
        Matrices.normalize(B);
        Matrices.normalize(pi);
    }
    public HMM (int N, int M, int maxIter) {
        this.N = N;
        this.M = M;
        this.maxIter = maxIter;
        A = new double[N][N];
        B = new double[N][M];
        pi = new double[N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = 1 / N + Math.random() * 0.01 / N;
            }
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                B[i][j] = 1 / N + Math.random() * 0.5 / N;
            }
        }
        for (int i = 0; i < N; i++) {
            pi[i] = 1 / N + Math.random() * 0.01 / N;
        }
        Matrices.normalize(A);
        Matrices.normalize(B);
        Matrices.normalize(pi);
    }



    public Integer[] getObs() {
        return obs;
    }

    public void load(InputStream stream, boolean haveObs) {
        /*
        Read three matrices (in this order);
        transition matrix, emission matrix, and initial state probability distribution.
        */
        Scanner scanner = new Scanner(stream);
        if (!scanner.hasNextInt()) {
            System.err.println("Can not read input file!");
            return;
        }
        int rowA, colA, rowB, colB, colPi;
        N = rowA = scanner.nextInt();
        colA = scanner.nextInt();
        A = new double[rowA][];
        Matrices.readMatrix(scanner, rowA, colA, A);
        rowB = scanner.nextInt();
        M = colB = scanner.nextInt();
        B = new double[rowB][];
        Matrices.readMatrix(scanner, rowB, colB, B);
        scanner.next();
        scanner.next();
        pi = new double[N];
        for (int i = 0; i < N; i++) {
            pi[i] = scanner.nextDouble();
        }
        if (haveObs) {
            T = scanner.nextInt();
            obs = new Integer[T];
            for (int i = 0; i < T; i++) {
                obs[i] = scanner.nextInt();
            }
        }

    }

    public double[] nextEmission(double[] stateDistribution) {
        /*
        Return the most likely emission at time(t)
         */
        if (A == null || B == null || stateDistribution == null) {
            System.out.println("Not initialized");
            return null;
        }
        // First find out after the transition how likely it is in state A, B, ...
        double[] afterTransition = Matrices.vecByMat(stateDistribution, A);

        // The find out how likely for each state to emit observation 1, 2, 3...
        double[] nextEmission = Matrices.vecByMat(afterTransition, B);
        return nextEmission;
    }

    public double[] estimateStateDistribution(Integer[] obs) {
        forwardAlgorithmEM(obs);
        return null;
    }

    public double forwardAlgorithm(Integer[] obs) {
        /*
        Estimate the probability of the made observation sequences by alpha-pass algorithm
         */
        // Initialize alpha vector
        double[] a = Matrices.vecByVec(pi, Matrices.getMatCol(B, obs[0]));
        double[][] alpha = new double[N][];
        for (int i = 0; i < N; i++) {
            alpha[i] = new double[obs.length];
        }

        Matrices.setMatCol(alpha, a, 0);

        /* Forward by calculating how likely to be in each state and then how likely to
         produce such emission in each state*/
        for (int i = 1; i < obs.length; i++) {
            a = Matrices.vecByMat(a, A);
            a = Matrices.vecByVec(a, Matrices.getMatCol(B, obs[i]));
            Matrices.setMatCol(alpha, a, i);
        }

        this.alpha = alpha;
        // Sum up all the likelihood to produce the observation in each state
        double sum = Arrays.stream(
                Matrices.getMatCol(alpha, obs.length-1)).sum();
        return sum;
    }

    public double[][] forwardAlgorithmEM(Integer[] obs) {
        /*
        Estimate the probability of the made observation sequences by alpha-pass algorithm
         */
        // Compute a0
        T = obs.length;
        double[][] alpha = new double[N][];
        Stack<Double> scales = new Stack<>();

        for (int i = 0; i < N; i++) {
            alpha[i] = new double[T];
        }
        Double scaleFactor = 0.;
        for (int i = 0; i < N; i++) {
            alpha[i][0] = pi[i] * B[i][obs[0]];
            scaleFactor += alpha[i][0];
        }
        scaleFactor = 1 / scaleFactor;
        scales.push(scaleFactor);
        // Scale the a0
        for (int i = 0; i < N; i++) {
            alpha[i][0] *= scaleFactor;
        }
        // Compute at(i)
        for (int t = 1; t < T; t++) {
            scaleFactor = 0.;
            for (int i = 0; i < N; i++) {
                alpha[i][t] = 0.;
                for (int j = 0; j < N; j++) {
                    alpha[i][t] += alpha[j][t-1] * A[j][i];
                }
                alpha[i][t] *= B[i][obs[t]];
                scaleFactor += alpha[i][t];
            }
            scaleFactor = 1 / scaleFactor;
            scales.push(scaleFactor);
            for (int i = 0; i < N; i++) {
                alpha[i][t] *= scaleFactor;
            }
        }
        this.alpha = alpha;
        this.scales = scales;
        return alpha;
    }

    public double[][] backwardAlgorithm(Integer[] obs) {
        // Initialization
        DecimalFormat df = new DecimalFormat("#.0000");
        double[][] beta = new double[N][];
        for (int i = 0; i < N; i++) {
            beta[i] = new double[obs.length];
        }
        for (int i = 0; i < N; i++) {
            beta[i][obs.length - 1] = 1.;
        }
        // Loop
        for (int i = obs.length - 2; i >= 0; i--) {
            for (int j = 0; j < N; j++) {
                Double prob = Matrices.elementWise(
                        Matrices.getMatRow(A, j),
                        Matrices.getMatCol(B, obs[i + 1]),
                        Matrices.getMatCol(beta, i+1));
                beta[j][i] = Double.valueOf(df.format(prob));
            }
        }
        this.beta = beta;
        Matrices.printMatrix(beta);
        return beta;
    }

    public void backwardAlgorithmEM(Integer[] obs) {
        // Initialization
        double[][] beta = new double[N][];
        double scaleFactor = scales.peek();

        for (int i = 0; i < N; i++) {
            beta[i] = new double[T];
        }
        for (int i = 0; i < N; i++) {
            beta[i][T - 1] = scaleFactor;
        }
        // The 尾-pass
        for (int t = T - 2; t >= 0; t--) {
            scaleFactor = scales.get(t);
            for (int i = 0; i < N; i++) {
                beta[i][t] = 0.;
                for (int j = 0; j < N; j++) {
                    beta[i][t] += A[i][j] * B[j][obs[t+1]] * beta[j][t+1];
                }
                beta[i][t] *= scaleFactor;
            }
        }
        this.beta = beta;
    }

    public void BWAlgorithm(Integer[] obs) {
        this.T = obs.length;
        int iter = 0;
        Double oldLogProb =  - Double.MAX_VALUE;
        // Initialization and memory allocation
        double[][][] gamma = new double[T][][];
        double[][] y = new double[T][];
        double[][] newA = new double[N][];
        double[][] newB = new double[N][];
        double[] newPi = new double[N];
        for (int i = 0; i < N; i++) {
            newA[i] = new double[N];
        }
        for (int i = 0; i < N; i++) {
            newB[i] = new double[M];
        }
        for (int i = 0; i < T; i++) {
            gamma[i] = new double[N][];
            for (int j = 0; j < N; j++) {
                gamma[i][j] = new double[N];
            }
        }
        for (int i = 0; i < T; i++) {
            y[i] = new double[N];
        }
        Double logProb;
        while (true) {
            forwardAlgorithmEM(obs);
            backwardAlgorithmEM(obs);
            // Compute di-gamma
            for (int t = 0; t < T - 1; t++) {
                for (int i = 0; i < N; i++) {
                    y[t][i] = 0.;
                    for (int j = 0; j < N; j++) {
                        Double product =  alpha[i][t] * A[i][j] * B[j][obs[t+1]] * beta[j][t+1];
                        y[t][i] += product;
                        gamma[t][i][j] = product;
                    }
                }
            }
            // Special case for y_{T-1}(i)
            for (int i = 0; i < N; i++) {
                y[T - 1][i] = alpha[i][T - 1];
            }

            /* Re-estimate A, B and pi */
            // Re-estimate pi
            for (int i = 0; i < N; i++) {
                newPi[i] = y[0][i];
            }

            // Re-estimate A
            for (int i = 0; i < N; i++) {
                Double denominator = 0.;
                // Denominator is the probability of being at state i
                for (int t = 0; t < T - 1; t++) {
                    denominator += y[t][i];
                }
                for (int j = 0; j < N; j++) {
                    Double numerator = 0.;
                    for (int t = 0; t < T - 1; t++) {
                        numerator += gamma[t][i][j];
                    }
                    newA[i][j] = numerator / denominator;
                }
            }

            // Re-estimate B
            for (int i = 0; i < N; i++) {
                Double denominator = 0.;
                for (int t = 0; t < T; t++) {
                    denominator += y[t][i];
                }
                for (int j = 0; j < M; j++) {
                    Double numerator = 0.;
                    for (int t = 0; t < T; t++) {
                        if (obs[t] == j) {
                            numerator += y[t][i];
                        }
                    }
                    newB[i][j] = numerator / denominator;
                }
            }

            // Compute logP(O|lambda)
            logProb = 0.;
            for (int t = 0; t < T; t++) {
                logProb += Math.log(scales.get(t));
            }
            logProb = - logProb;
            A = newA.clone();
            B = newB.clone();
            pi = newPi.clone();

            iter += 1;
            if (Math.abs(logProb - oldLogProb) < 1E-300 || iter > maxIter) {
//                System.out.println("Iterations end at:" + iter);
//                System.out.println("Found!Output params:");
//                System.out.println("A:");
//                Matrices.printMatrix(A);
//                System.out.println("\nB:");
//                Matrices.printMatrix(B);
//                System.err.println("\npi:");
//                Matrices.printVector(pi);
//                System.out.println("what1");
                break;
            }  else {
                oldLogProb = logProb;
            }
        }

    }
    class Op {
        double numerator = 0.;
        double denominator = 0.;
        void add(double v1, double v2) {
            numerator += v1;
            denominator += v2;
        }
        void setZero() {
            numerator = 0.;
            denominator = 0.;
        }
    }
    public void trainMultiSeq(ArrayList<Integer[]> seq) {
        double[] specialPi = new double[N];
        Op[][] opA = new Op[N][];
        for (int i = 0; i < N; i++) {
            opA[i] = new Op[N];
            for (int j =0; j < N; j++) {
                opA[i][j] = new Op();
            }
        }
        Op[][] opB = new Op[N][];
        for (int i = 0; i < N; i++) {
            opB[i] = new Op[M];
            for (int j = 0; j < M; j++) {
                opB[i][j] = new Op();
            }
        }
        int iter = 0;
        double oldLogProb =  - Double.MAX_VALUE;
        while (true) {
            double overallLogProb = 0.;
            // 浠ヤ笅for寰幆瀹屾垚涓�娆¤缁�
            for (int l = 0; l < seq.size(); l++) {
                Integer[] obs = seq.get(l);
                int T = obs.length;
                double[][][] gamma = new double[T][][];
                double[][] y = new double[T][];
                for (int i = 0; i < T; i++) {
                    gamma[i] = new double[N][];
                    for (int j = 0; j < N; j++) {
                        gamma[i][j] = new double[N];
                    }
                }
                for (int i = 0; i < T; i++) {
                    y[i] = new double[N];
                }
                forwardAlgorithmEM(obs);
                backwardAlgorithmEM(obs);
                // Compute di-gamma
                for (int t = 0; t < T - 1; t++) {
                    for (int i = 0; i < N; i++) {
                        y[t][i] = 0.;
                        for (int j = 0; j < N; j++) {
                            Double product =  alpha[i][t] * A[i][j] * B[j][obs[t+1]] * beta[j][t+1];
                            y[t][i] += product;
                            gamma[t][i][j] = product;
                        }
                    }
                }
                // Special case for y_{T-1}(i)
                for (int i = 0; i < N; i++) {
                    y[T - 1][i] = alpha[i][T - 1];
                }
                // Re-estimate pi
                for (int i = 0; i < N; i++) {
                    specialPi[i] += y[0][i];
                }

                // Re-estimate A
                for (int i = 0; i < N; i++) {
                    double denominator = 0.;
                    // Denominator is the probability of being at state i
                    for (int t = 0; t < T - 1; t++) {
                        denominator += y[t][i];
                    }
                    for (int j = 0; j < N; j++) {
                        double numerator = 0.;
                        for (int t = 0; t < T - 1; t++) {
                            numerator += gamma[t][i][j];
                        }
                        opA[i][j].add(numerator, denominator);
                    }
                }
                // Re-estimate B
                for (int i = 0; i < N; i++) {
                    double denominator = 0.;
                    for (int t = 0; t < T; t++) {
                        denominator += y[t][i];
                    }
                    for (int j = 0; j < M; j++) {
                        double numerator = 0.;
                        for (int t = 0; t < T; t++) {
                            if (obs[t] == j) {
                                numerator += y[t][i];
                            }
                        }
                        opB[i][j].add(numerator, denominator);
                    }
                }
                // Compute logP(O|lambda)
                double logProb = 0.;
                for (int t = 0; t < T; t++) {
                    logProb += Math.log(scales.get(t));
                }
                logProb = - logProb;
                overallLogProb += logProb;
            } // 缁撴潫涓�娆¤缁�
            // 姹備笁涓弬鏁癆 B pi鐨勫钩鍧囧��
            for (int i = 0; i < N; i++) {
                pi[i] = specialPi[i] / seq.size();
            }
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    A[i][j] = opA[i][j].numerator / opA[i][j].denominator;
                }
            }
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < M; j++) {
                    B[i][j] = opB[i][j].numerator / opB[i][j].denominator;
                }
            }
            // 姹傚畬鍊间箣鍚庯紝閲嶇疆绱姞鍣�
            for (int i = 0; i < N; i++) {
                specialPi[i] = 0.;
            }
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    opA[i][j].setZero();
                }
                for (int j = 0; j < M; j++) {
                    opB[i][j].setZero();
                }
            }
            // 鍒ゆ柇鏀舵暃鏉′欢
            iter += 1;
            if (iter > maxIter || Math.abs(overallLogProb - oldLogProb) < 1E-300) {
                break;
            } else {
                oldLogProb = overallLogProb;
            }
        }
    }
    
    public void estimateStateSeq(Integer[] obs) {
        /*
        Implement Viterbi algorithm to compute the most likely sequence of hidden states given
        the observations.
         */
        // Initialization
        // Record the optimal previous state for the current state in each time step
        Integer[][] preState = new Integer[N][];

        // Record the probability of the corresponding state in preState
        double[][] T = new double[N][];

        // Allocate memory for each row in the matrix
        for (int i = 0; i < N; i++) {
            preState[i] = new Integer[obs.length];
            T[i] = new double[obs.length];
        }

        // Intialize by multiplying pi and column of B
        for (int i = 0; i < N; i++) {
            preState[i][0] = 0;
            T[i][0] = B[i][obs[0]] * pi[i];
        }

        // Delta holds the result column in each loop
        double[] delta = Matrices.vecByVec(pi, Matrices.getMatCol(B, obs[0]));

        /* Whatever state u are in(Assuming u are in A), consider the likelihood that
             you just went from B, C, D...then transfer to A and emit the given result;
             then consider you are in B and do the same thing, then C, D...*/
        double[] newProb, newDelta = new double[N];
        for (int i = 1; i < obs.length; i++) {
            for (int j = 0; j < N; j++) {
                double[] emitProb = new double[N];
                Arrays.fill(emitProb, B[j][obs[i]]);
                newProb = Matrices.vecByVec(Matrices.vecByVec(delta, Matrices.getMatCol(A, j)), emitProb);
                double maxProb = Arrays.stream(newProb).max().getAsDouble();
                // 杩欓噷鏈塨ug锛岄渶瑕佹敼
                preState[j][i] = Arrays.asList(newProb).indexOf(maxProb);
                newDelta[j] = maxProb;
                T[j][i] = maxProb;
            }
            delta = newDelta.clone();
        }

        // Compute the path by tracing back from the end
        LinkedList<Integer> path = new LinkedList<>();
        Integer end;
        double[] lastCol = Matrices.getMatCol(T, obs.length - 1);
        int terminal = Arrays.asList(lastCol).indexOf(Arrays.stream(lastCol).max().getAsDouble());
        path.push(terminal);
        for (int i = obs.length - 1; i > 0; i--) {
            path.push(preState[path.peekFirst()][i]);
        }
        Iterator iterator = path.iterator();
        while (iterator.hasNext()) {
            System.out.print(iterator.next() + " ");
        }
        System.out.println();
    }

//    public static void main(String[] args) {
//        double[][] m = new double[][]{{1, 2, 3}, {4, 5, 6}};
//        Matrices.printMatrix(m);
//        Matrices.normalize(m);
//        System.err.println();
//        Matrices.printMatrix(m);
//    }
}
