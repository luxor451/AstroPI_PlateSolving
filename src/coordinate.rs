pub struct Arcdegrees {
    pub degrees: i64,    // in degrees
    pub arcminutes: i64, // in minutes
    pub arcseconds: f64, // in seconds
}

pub struct RaHoursMinutesSeconds {
    pub hours: i64,   // in hours
    pub minutes: i64, // in minutes
    pub seconds: f64, // in seconds
}

pub struct CoordinateEquatorial {
    pub ra: RaHoursMinutesSeconds, // in degrees
    pub dec: Arcdegrees,           // in degrees
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

    #[allow(warnings)]
    pub fn from_radians(radians: f64) -> Self {
        let degrees_total = radians.to_degrees();
        let degrees = degrees_total.trunc() as i64;
        let minutes_total = (degrees_total - degrees as f64) * 60.0;
        let arcminutes = minutes_total.trunc() as i64;
        let arcseconds = (minutes_total - arcminutes as f64) * 60.0;
        Arcdegrees {
            degrees,
            arcminutes,
            arcseconds,
        }
    }
}

impl std::fmt::Display for Arcdegrees {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}° {}' {:.2}\"",
            self.degrees, self.arcminutes, self.arcseconds
        )
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
    #[allow(warnings)]
    pub fn to_arcdegrees(&self) -> Arcdegrees {
        let total_seconds = (self.hours * 3600 + self.minutes * 60) as f64 + self.seconds * 60.0;
        let degrees = total_seconds / 240.0; // 1 hour = 15 degrees, so 1 second = 15/3600 degrees
        let arcminutes = (degrees.fract() * 60.0) as i64;
        let arcseconds = (degrees.fract() * 3600.0).fract() * 60.0;

        Arcdegrees::new(degrees as i64, arcminutes, arcseconds)
    }

    pub fn to_degrees(&self) -> f64 {
        (360.0 / 24.0)
            * (self.hours as f64 + (self.minutes as f64 / 60.0) + (self.seconds / 3600.0))
    }

    pub fn to_hours(&self) -> f64 {
        self.hours as f64 + (self.minutes as f64 / 60.0) + (self.seconds / 3600.0)
    }

    pub fn to_radians(&self) -> f64 {
        self.to_degrees() * std::f64::consts::PI / 180.0
    }
    #[allow(warnings)]
    pub fn from_radians(radians: f64) -> Self {
        let degrees = radians.to_degrees();
        let total_hours = degrees / 15.0;
        let hours = total_hours.trunc() as i64;
        let total_minutes = (total_hours - hours as f64) * 60.0;
        let minutes = total_minutes.trunc() as i64;
        let seconds = (total_minutes - minutes as f64) * 60.0;
        RaHoursMinutesSeconds {
            hours,
            minutes,
            seconds,
        }
    }
}

impl std::fmt::Display for RaHoursMinutesSeconds {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}h {}m {:.2}s", self.hours, self.minutes, self.seconds)
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
    #[allow(warnings)]
    pub fn to_radians(&self) -> (f64, f64) {
        let ra_rad = self.ra.to_radians();
        let dec_rad = self.dec.to_radians();
        (ra_rad, dec_rad)
    }

    pub fn from_radians(ra_rad: f64, dec_rad: f64) -> Self {
        // Convert RA from radians to total hours
        let total_ra_hours = ra_rad * 12.0 / std::f64::consts::PI;
        let ra_h = total_ra_hours.trunc() as i64;
        let ra_m_total = total_ra_hours.fract() * 60.0;
        let ra_m = ra_m_total.trunc() as i64;
        let ra_s = ra_m_total.fract() * 60.0;
        let ra = RaHoursMinutesSeconds::new(ra_h, ra_m, ra_s);

        // Convert Dec from radians to total degrees
        let total_dec_degrees = dec_rad.to_degrees();
        let dec_d = total_dec_degrees.trunc() as i64;
        let dec_m_total = total_dec_degrees.abs().fract() * 60.0;
        let dec_m = dec_m_total.trunc() as i64;
        let dec_s = dec_m_total.fract() * 60.0;
        let dec = Arcdegrees::new(dec_d, dec_m, dec_s);

        CoordinateEquatorial { ra, dec }
    }

    /// Creates a new `CoordinateEquatorial` from RA and Dec in degrees.
    /// 
    /// # Arguments
    /// 
    /// * `ra_deg` - Right Ascension in degrees (0-360).
    /// * `dec_deg` - Declination in degrees (-90 to +90).
    pub fn from_degrees(ra_deg: f64, dec_deg: f64) -> Self {
        // Convert RA from degrees to hours (15 degrees = 1 hour)
        let total_ra_hours = ra_deg / 15.0;
        let ra_h = total_ra_hours.trunc() as i64;
        let ra_m_total = total_ra_hours.fract().abs() * 60.0;
        let ra_m = ra_m_total.trunc() as i64;
        let ra_s = ra_m_total.fract() * 60.0;
        let ra = RaHoursMinutesSeconds::new(ra_h, ra_m, ra_s);

        // Convert Dec from degrees
        let dec_d = dec_deg.trunc() as i64;
        let dec_m_total = dec_deg.abs().fract() * 60.0;
        let dec_m = dec_m_total.trunc() as i64;
        let dec_s = dec_m_total.fract() * 60.0;
        let dec = Arcdegrees::new(dec_d, dec_m, dec_s);

        CoordinateEquatorial { ra, dec }
    }
}

impl std::fmt::Display for CoordinateEquatorial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "RA: {}, Dec: {}", self.ra, self.dec)
    }
}
