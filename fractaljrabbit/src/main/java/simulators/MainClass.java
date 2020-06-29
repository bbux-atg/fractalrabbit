package simulators;

/** THIS PUBLIC DOMAIN SOFTWARE WAS PRODUCED BY AN EMPLOYEE OF U.S. GOVERNMENT
 * Produces an output csv file where each line has the format
 * identifier (integer), time (days), latitude, longitude
 * or else
 * identifier (integer), time (days), x (km), y (km)
 * Here identifier refers to a specific traveler.
 *
 * A set of input parameters is here given explicitly.
 * Later these will be input through a csv file
 * args[0] will be a path to write the output file
 *
 * Extra functionality is the generation of MULTIPLE trajectories in the same
 * set of points, with some co-travelers
 *
 * Ran successfully Februay 1, 2019
 *
 * April 3, 2019: modified to read in 7 principal parameters via parameters.csv input file
 */

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.csv.CSVRecord;
import utilities.EuclideanPoint;
import utilities.PoissonVariate;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author RWRD
 *
 */
public class MainClass {
    /*
     * Three tiers of the FRACTALRABBIT simulator
     */
    private final AgoraphobicPoints app;
    private final Retropreferential rpp;
    private final SporadicReporter spore;
    private final List<List<Integer>> trajectoryList;
    private final Map<Integer, Integer> trajectoryAssignment; // traveler to trajectory number
    private final List<Integer> coTravellers;
    private final Map<Integer, List<Double>> reportTimesAssignment; // traveler to list of report times
    private final Map<Integer, List<Integer>> reportPlacesAssignment; // traveler to list of point references (at given times)

    /*
     * Minimal parameters to create an instance of each simulator Some are
     * hard-coded, but seven are in the parameters.csv file
     */
    public MainClass(Parameters params) {
        this.app = new AgoraphobicPoints(params.n, params.dimension, params.h, params.theta);
        EuclideanPoint[] euclidPoints = new EuclideanPoint[app.getPoints().size()];
        Arrays.setAll(euclidPoints, j -> new EuclideanPoint(app.getPoints().get(j)));
        this.rpp = new Retropreferential(euclidPoints, params.exponent);
        this.spore = new SporadicReporter(euclidPoints, params.speedbound, params.kilometersPerUnit);
        this.trajectoryList = new ArrayList<>();
        this.coTravellers = new ArrayList<>();
        /*
         * The following three maps all have the SAME key set, i.e. shuffled travelers
         */
        this.trajectoryAssignment = new HashMap<>(); // assigns traveler to trajectory
        this.reportTimesAssignment = new HashMap<>(); // assigns traveler to times
        this.reportPlacesAssignment = new HashMap<>(); // assigns traveler to places
    }

    /*
     * Assign travelers to trajectories, with randomization
     */
    public void travelerConfigure(int numCoTravelers, int numTravelers) {
        Integer[] travelerIDs = new Integer[numTravelers];
        Arrays.setAll(travelerIDs, i -> i);
        /*
         * Shuffle the traveler IDs, and assign the first #(numCoTravelers) to first
         * trajectory
         */
        List<Integer> travelerShuffle = Arrays.asList(travelerIDs);
        Collections.shuffle(travelerShuffle);
        // All the co-travelers are assigned to trajectory 0, and their labels are
        // recorded
        for (int t = 0; t < numCoTravelers; t++) {
            this.trajectoryAssignment.put(travelerShuffle.get(t), 0);
            this.coTravellers.add(travelerShuffle.get(t));
        }
        // Other travelers are assigned to trajectories 1 through numTrajectories - 1
        for (int t = numCoTravelers; t < numTravelers; t++) {
            this.trajectoryAssignment.put(travelerShuffle.get(t), t - numCoTravelers + 1);
        }
    }

    public static void main(String[] argv) throws IOException {

        Args args = new Args();
        JCommander.newBuilder().addObject(args).build().parse(argv);

        String csvOutputFile = args.outFile;

        Parameters params;
        if (args.csvParams.isBlank()) {
            params = fromArgs(args);
        } else {
            try {
                params = fromCsv(args.csvParams);
            } catch (IOException e) {
                System.out.println("Unable to load params from: " + args.csvParams);
                return;
            }
        }

        /*
         * Instantiate generator with these parameters. Here the points are generated.
         */
        MainClass frg = new MainClass(params);
        System.out.println(params.n + " points generated for fractal dimension " + params.h);
        /*
         * Generate the trajectories (points exist already). Each trajectory needs a
         * random phi, random # of steps, random start point,
         */
        double phiMean = 10.0; // Retro-preferential exploration - mean value
        PoissonVariate phiRandom = new PoissonVariate(phiMean);
        double stepsMean = 100.0; // Retro-preferential trajectory length - mean value
        PoissonVariate stepsRandom = new PoissonVariate(stepsMean);
        Random g = new Random();
        /*
         * Set up co-travel structure
         */
        frg.travelerConfigure(params.numCoTravelers, params.numTravelers);
        int numTrajectories = params.numTravelers - params.numCoTravelers + 1;// all co-travelers use same trajectory
        // Loop through trajectories
        for (int t = 0; t < numTrajectories; t++) {
            // new trajectory with Poisson phi, Poisson # of steps, uniform random start
            // point
            frg.trajectoryList.add(frg.rpp.trajectory(phiRandom.generate(), stepsRandom.generate(), g.nextInt(params.n)));
        }

        /*
         * Sporadic reporting parameters
         */
        double delta = 0.001; // Sporadic Reporter - Pareto cutoff
        double alpha = -1.5; // Sporadic Reporter - Pareto tail
        PoissonVariate countRandom = new PoissonVariate(params.countMean);

        System.out.println(numTrajectories + " trajectories generated, in which travelers assigned to same trajectory are: ");
        for (Integer p : frg.coTravellers) {
            System.out.print(p + ", ");
        }
        System.out.println();
        /*
         * Loop through travelers, attaching trajectory to each. p < # travelers,
         * defined above as (numTrajectories - 1 + numCoTravelers) Co-travelers will
         * receive same trajectory.
         *
         * ERROR HERE!
         */
        int currentTrajectory;
        for (Integer p : frg.trajectoryAssignment.keySet()) {
            currentTrajectory = frg.trajectoryAssignment.get(p); // trajectory # for traveler p
            frg.spore.embedTrajectory(frg.trajectoryList.get(currentTrajectory)); // embed in continuous time
            frg.spore.generateReports(frg.trajectoryList.get(currentTrajectory), countRandom.generate(), delta, alpha, params.days); // generate sporadic reports for this trajectory
            frg.reportTimesAssignment.put(p, List.copyOf(frg.spore.getReportTimes())); // tag p with Report Times
            frg.reportPlacesAssignment.put(p, List.copyOf(frg.spore.getReportPlaces())); // tag p with Report Places
        }
        /*
         * Check that the trajectories differ , by inspecting first element of each -
         * FAILED!
         */
        for (Integer p : frg.reportTimesAssignment.keySet()) {
            System.out
                .println(
                    "First report for traveler "
                        + p
                        + " at time "
                        + frg.reportTimesAssignment.get(p).get(0)
                        + " is place "
                        + frg.reportPlacesAssignment.get(p).get(0));
        }
        System.out.println("Reports have been generated for " + params.numTravelers + " travelers.");
        /*
         * Need only format the data as (traveler ID, time, x in km, y in km), and write
         * to csv using Apache Commons CSV
         */
        try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(csvOutputFile + ".csv"));
            CSVPrinter csvPrinter = new CSVPrinter(writer, CSVFormat.DEFAULT.withHeader("ID", "Days", "x(km)", "y(km)"))) {
            double x, y, time;
            int place;
            int lines = 0;
            // Write to file ERROR - repetitions of one report set!

            for (Integer p : frg.reportTimesAssignment.keySet()) {
                for (int t = 0; t < frg.reportTimesAssignment.get(p).size(); t++) {
                    time = frg.reportTimesAssignment.get(p).get(t);
                    place = frg.reportPlacesAssignment.get(p).get(t);
                    if (t == 0) {
                        System.out
                            .println(
                                "First report for traveler "
                                    + p
                                    + " is point "
                                    + frg.reportPlacesAssignment.get(p).get(t)
                                    + " at time "
                                    + frg.reportTimesAssignment.get(p).get(t));
                    }
                    x = params.kilometersPerUnit * frg.app.getPoints().get(place)[0];
                    y = params.kilometersPerUnit * frg.app.getPoints().get(place)[1];
                    csvPrinter.printRecord(p, time, x, y);
                    lines++;
                }
            }

            System.out.println("CSV file created with " + lines + " lines, called " + csvOutputFile + ".csv");
            csvPrinter.flush();
            writer.flush();
        }

        /*
         * SECOND VERSION! (traveler ID, time, place ID), and write to csv using Apache
         * Commons CSV
         */
        try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(csvOutputFile + "PLACES.csv"));
            CSVPrinter csvPrinter = new CSVPrinter(writer, CSVFormat.DEFAULT.withHeader("Traveler ID", "Days", "Place ID"))) {
            int lines = 0;
            // Write to file ERROR - repetitions of one report set!

            for (Integer p : frg.reportTimesAssignment.keySet()) {
                for (int t = 0; t < frg.reportTimesAssignment.get(p).size(); t++) {
                    csvPrinter.printRecord(p, frg.reportTimesAssignment.get(p).get(t), frg.reportPlacesAssignment.get(p).get(t));
                    lines++;
                }
            }
            System.out.println("CSV Places file created with " + lines + " lines, called " + csvOutputFile + "PLACES.csv");
            csvPrinter.flush();
            writer.flush();
        }
    }

    private static Parameters fromCsv(String inputFile) throws IOException {
        List<String> par = new ArrayList<>(); //parameters as strings
        try(Reader inputParameterFile = new FileReader(inputFile)) {
            Iterable<CSVRecord> records = CSVFormat.RFC4180.parse(inputParameterFile);
            for (CSVRecord record : records) {
                System.out.println(record.get(0) + " has value " + record.get(1));
                par.add(record.get(1).replaceAll("\\s", "")); // remove whitespace
            }
        }

        Parameters params = new Parameters();
        params.dimension = Integer.parseInt(par.get(0));
        params.h = Double.parseDouble(par.get(1));
        params.n = Integer.parseInt(par.get(2));
        params.numTravelers = Integer.parseInt(par.get(3));
        params.numCoTravelers = Integer.parseInt(par.get(4));
        params.days = Double.parseDouble(par.get(5));
        params.countMean = Double.parseDouble(par.get(6));

        return params;
    }

    private static Parameters fromArgs(Args args) {
        Parameters params = new Parameters();
        params.dimension = args.spaceDemensions;
        params.h = args.fractalDimension;
        params.n = args.numPoints;
        params.numTravelers = args.numTravellers;
        params.numCoTravelers = args.numCoTravellers;
        params.days = args.numDays;
        params.countMean = args.numReports;

        return params;
    }

    private static class Parameters {
        /*
         * Small number of parameters which will often change
         */
        int dimension;
        double h; // fractal dimension for AgoraphobicPoints
        int n;// # points in Agoraphobic points process
        int numTravelers;
        int numCoTravelers;
        double days; // Sporadic Reporter - duration in days
        double countMean; // Sporadic Reporter - # reports - mean value

        /*
         * Large number of parameters which seldom change
         */
        double theta = 0.75; // restart rate for AgoraphobicPoints
        double exponent = -2.0; // convert distances to Retro-preferential matrix
        double speedbound = 50.0; // Sporadic Reporter - units per day
        double kilometersPerUnit = 200.0; // Sporadic Reporter - # km per one distance unit for points
    }

    private static class Args {
        @Parameter(names = { "-o", "--output" }, description = "Prefix for output file names", required = true)
        private String outFile;

        @Parameter(names = { "-p", "--csv-params" }, description = "Legacy CSV params file")
        private String csvParams = "";

        @Parameter(names = { "-sd", "--space-dimensions" }, description = "Number of Space Dimensions")
        private Integer spaceDemensions = 2;

        @Parameter(names = { "-fd", "--fractal-dimensions" }, description = "fractal dimension for AgoraphobicPoints")
        private Double fractalDimension = 1.33;

        @Parameter(names = { "-np", "--num-points" }, description = "Number of points in Agoraphobic points process")
        private Integer numPoints = 2;

        @Parameter(names = { "-nt", "--num-travellers" }, description = "Number of travellers")
        private Integer numTravellers = 1000;

        @Parameter(names = { "-nc", "--num-co-travellers" }, description = "Number of co-travellers")
        private Integer numCoTravellers = 2;

        @Parameter(names = { "-nd", "--num-days" }, description = "Number of days")
        private Double numDays = 60.0;

        @Parameter(names = { "-nr", "--num-reports" }, description = "Number of reports")
        private Double numReports = 200.0;
    }
}
