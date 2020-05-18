
   Volume<1> volume( const Point<1>& p0, const Point<1>& p1 )
  {
      const Direction<1> dir=p1-p0;
      Volume<1> vol;
      vol.volume=length(dir);
      vol.centre=p0+0.5*dir;
      return vol;
  }
