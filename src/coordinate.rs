pub struct Arcdegrees {
    pub degrees: i64, // in degrees
    pub arcminutes: i64, // in minutes
    pub arcseconds: f64, // in seconds
}

pub struct RaHoursMinutesSeconds {
    pub hours: i64, // in hours
    pub minutes: i64, // in minutes
    pub seconds: f64, // in seconds
}

pub struct CoordinateEquatorial {
    pub ra: RaHoursMinutesSeconds, // in degrees
    pub dec: Arcdegrees, // in degrees
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Star {
    /// Right Ascension in degrees.
    pub ra: f64,
    /// Declination in degrees.
    pub dec: f64,
}

impl Arcdegrees {
    pub fn new(degrees: i64, arcminutes: i64, arcseconds: f64) -> Self {
        Arcdegrees {
            degrees,
            arcminutes,
            arcseconds,
        }
    }

    pub fn to_degrees(&self) -> f64 {
        self.degrees as f64 + self.arcminutes as f64 / 60.0 + self.arcseconds / 3600.0
    }

    pub fn to_radians(&self) -> f64 {
        self.to_degrees() * std::f64::consts::PI / 180.0
    }
}

impl RaHoursMinutesSeconds {
    pub fn new(hours: i64, minutes: i64, seconds: f64) -> Self {
        RaHoursMinutesSeconds {
            hours,
            minutes,
            seconds,
        }
    }

    pub fn to_arcdegrees(&self) -> Arcdegrees {
        let total_seconds = (self.hours * 3600 + self.minutes * 60 + self.seconds as i64) as f64;
        let degrees = total_seconds / 240.0; // 1 hour = 15 degrees, so 1 second = 15/3600 degrees
        let arcminutes = (degrees.fract() * 60.0) as i64;
        let arcseconds = (degrees.fract() * 3600.0).fract() * 60.0;

        Arcdegrees::new(degrees as i64, arcminutes, arcseconds)
    }

    pub fn to_degrees(&self) -> f64 {
        let total_seconds = (self.hours * 3600 + self.minutes * 60) as f64 + self.seconds * 60.0;
        total_seconds / 240.0 // 1 hour = 15 degrees, so 1 second = 15/3600 degrees
    }

    pub fn to_radians(&self) -> f64 {
        self.to_degrees() * std::f64::consts::PI / 180.0
    }
}

impl CoordinateEquatorial {
    pub fn new(ra: RaHoursMinutesSeconds, dec: Arcdegrees) -> Self {
        CoordinateEquatorial { ra, dec }
    }

    pub fn to_degrees(&self) -> (f64, f64) {
        let ra_deg = self.ra.to_degrees();
        let dec_deg = self.dec.to_degrees();
        (ra_deg, dec_deg)
    }

    pub fn to_radians(&self) -> (f64, f64) {
        let ra_rad = self.ra.to_radians();
        let dec_rad = self.dec.to_radians();
        (ra_rad, dec_rad)
    }
}