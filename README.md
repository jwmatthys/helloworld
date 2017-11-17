# swe~ - Swiss Ephemeris in Pd #

Swiss Ephemeris is a function package of astronomical calculations that serves the needs of astrologers, archaeoastronomers, and, depending on purpose, also the needs of astronomers. It includes long-term ephemerides for the Sun, the Moon, the planets, more than 300’000 asteroids, historically relevant fixed stars and some “hypothetical” objects.

ftp://ftp.astro.ch/pub/swisseph/

The core part of Swiss Ephemeris is a compression of the JPL-Ephemeris DE431, which covers roughly the time range 13000 BCE to 17000 CE. Additional files availablle from ftp://ftp.astro.ch/pub/swisseph/ephe can be downloaded to expand the time range and/or access additional astronomical objects (include 300,000+ asteroids).

## Use ##

### Astronomical Body ###

swe~ can be initialized with a single argument, a number representing a planet or astronomical body.

    0 = Sun
    1 = Moon
    2 = Mercury
    3 = Venus
    4 = Mars
    5 = Jupiter
    6 = Saturn
    7 = Uranus
    8 = Neptune
    9 = Pluto
 
Additional bodies are listed at ftp://ftp.astro.ch/pub/swisseph/doc/swephprg.htm#_Toc471829059

The body can also be set with a ```p``` message:

    p 5

### Dates ###

In general swe~ uses dates in Julian calendar format. You can set the date with a ```b``` message (for Gregorian dates) or ```j``` for Julian dates:

    b 1974 6 27 0.25    (year, month, day, fraction of day)
    
    j 2245364.5
    
For best results, enclose Julian dates in quotes (single or double) to get full double accuracy:

    j '2442225'

### Getting output ###

Once the date and body are set, send a bang to the inlet to get the results.

From left to right, the outlets provide:

    [signal] [longitude] [latitude] [distance] [Julian yr] [Error bang]

The Julian year is output as a symbol to maintain accuracy. It can be converted back into a float value with a [float] object.

### Stepping forward in time ###

Send a float value to the inlet to increment time (in Julian days). The new data is output after each float.

The ```step``` message is used for automatic time increments used with audio signals, loops, and arrays. (See below.)

    step 0.5

### IFLAGS ###

swe~ has many additional options beyond the defaults. These can be set with an ```iflag``` message. For instance:

    iflag HELCTR (heliocentric position)
    iflag RADIANS (return results in radians rather than degrees)
    iflag DEFAULT (reset to default options)

The various IFLAGS are described here: ftp://ftp.astro.ch/pub/swisseph/doc/swephprg.htm#_Toc471829060

### Topocentric position ###

By default swe~ returns values calculated from the center of the earth. For topocentric position, send a ```topo``` message followed by latitude (degrees), longitude (degrees), and altitude (meters):

    topo 43.044 -88.23 252
    
The ```topo``` message resets all IFLAGS to defaults and activates the TOPOCTR flag. Send iflag DEFAULT to reset.

### Path ###

Set the path to additional ephemeris files with ```path```. Absolute paths only, no relative or abbreviated paths.

### Automatically filling arrays ###

swe~ can automatically fill Pd arrays with incremental data with the ```array``` message.

    array myarrayname "2442225" "2443225" long  (array name, start date, end date, output type)
    
The array is automatically resized to contain all the data, based on the time span and step size.

Options for output type are:

    long - longitude (default)
    lat - latitude
    dist - distance
    sin - sin(longitude in radians) (good for audio, since values are in range -1 to 1)

Be sure to set the body and step before attempting to fill the array.

## Calculating signals ##

_swe~ can be CPU-intensive. Often it is not feasible to calculate audio signals in real time._

By default the signal outlet of swe~ is deactivated. Upon receiving an ```audio 1``` message, swe~ will begin incrementing the date by the step size. The signal output is the sin of the longitude in radians, so it never exceeds the range -1 to 1.

Since swe~ is calculating so many values each second, it will usually exceed the built-in ephemeris quickly. Once the data is out of range, a message appears in the console and the signal outlet is deactivated. You can download additional ephemeris files at ftp://ftp.astro.ch/pub/swisseph/ephe to extend the working range.

### Looping ###

Since continuous incremental calculation at audio rate rapidly outruns the edges of the ephemeris, it is often better to loop a time range.

    loop "2442225" "2443225" 10  (start date, end date, step)
    
After setting the loop parameters, you must activate the audio with ```audio 1``` to begin the signal.

Values set for the loop do not overwrite the date and step size set for non-audio rate calculations.
