
namespace geom
{
   template<floating_point Real>
   Volume<1,Real> volume( const Point<1,Real>& p0, const Point<1,Real>& p1 )
  {
      const Direction<1,Real> dir=p1-p0;
      Volume<1,Real> vol;
      vol.volume=length(dir);
      vol.centre=p0+0.5*dir;
      return vol;
  }
}
