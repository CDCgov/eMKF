import os
import json
import requests
from datetime import datetime
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd

# Configuration
GITHUB_TOKEN = os.environ['GITHUB_TOKEN']
REPOSITORY = os.environ['REPOSITORY']
TRAFFIC_DIR = Path('.traffic')
HEADERS = {
    'Authorization': f'token {GITHUB_TOKEN}',
    'Accept': 'application/vnd.github.v3+json'
}

def fetch_traffic_data(endpoint):
    """Fetch traffic data from GitHub API"""
    url = f'https://api.github.com/repos/{REPOSITORY}/traffic/{endpoint}'
    response = requests.get(url, headers=HEADERS)
    response.raise_for_status()
    return response.json()

def fetch_stargazers():
    """Fetch all stargazers with timestamps"""
    url = f'https://api.github.com/repos/{REPOSITORY}/stargazers'
    headers = {
        'Authorization': f'token {GITHUB_TOKEN}',
        'Accept': 'application/vnd.github.v3.star+json'
    }
    
    all_stars = []
    page = 1
    per_page = 100
    
    while True:
        response = requests.get(
            url,
            headers=headers,
            params={'page': page, 'per_page': per_page}
        )
        response.raise_for_status()
        
        stars = response.json()
        if not stars:
            break
        
        all_stars.extend(stars)
        page += 1
        
        # Stop if we got fewer results than requested (last page)
        if len(stars) < per_page:
            break
    
    return all_stars

def aggregate_stars_by_date(stargazers):
    """Aggregate stargazers by date"""
    daily_stars = defaultdict(int)
    
    for star in stargazers:
        date = star['starred_at'][:10]  # Get YYYY-MM-DD
        daily_stars[date] += 1
    
    # Convert to cumulative count and format for storage
    sorted_dates = sorted(daily_stars.keys())
    cumulative = 0
    result = []
    
    for date in sorted_dates:
        cumulative += daily_stars[date]
        result.append({
            'timestamp': date,
            'count': daily_stars[date],  # Daily new stars
            'total': cumulative  # Total stars up to this date
        })
    
    return result

def load_historical_data(filename):
    """Load historical data from JSON file"""
    filepath = TRAFFIC_DIR / filename
    if filepath.exists():
        with open(filepath, 'r') as f:
            return json.load(f)
    return []

def merge_traffic_data(historical, new_data):
    """Merge new traffic data with historical data, avoiding duplicates"""
    existing_dates = {item['timestamp'] for item in historical}
    
    for item in new_data:
        if item['timestamp'] not in existing_dates:
            historical.append(item)
    
    # Sort by timestamp
    historical.sort(key=lambda x: x['timestamp'])
    return historical

def save_data(filename, data):
    """Save data to JSON file"""
    TRAFFIC_DIR.mkdir(exist_ok=True)
    filepath = TRAFFIC_DIR / filename
    with open(filepath, 'w') as f:
        json.dump(data, indent=2, fp=f)

def create_ascii_chart(data, metric_name, max_width=60):
    """Create ASCII chart from traffic data"""
    if not data:
        return f"No {metric_name} data available\n"
    
    # Get last 30 data points
    recent_data = data[-30:]
    
    dates = [item['timestamp'][:10] for item in recent_data]
    values = [item.get('count', 0) for item in recent_data]
    
    if not values or max(values) == 0:
        return f"No {metric_name} data available\n"
    
    max_val = max(values)
    chart = f"\n{metric_name.upper()} (Last 30 days)\n"
    chart += "=" * (max_width + 20) + "\n"
    
    for date, val in zip(dates, values):
        bar_len = int((val / max_val) * max_width) if max_val > 0 else 0
        bar = '█' * bar_len
        chart += f"{date} | {bar} {val}\n"
    
    chart += "=" * (max_width + 20) + "\n"
    return chart

def create_plot(views_data, clones_data, stars_data):
    """Create matplotlib plot of traffic data"""
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))
    
    # Views plot
    if views_data:
        df_views = pd.DataFrame(views_data)
        df_views['date'] = pd.to_datetime(df_views['timestamp'])
        
        ax1.plot(df_views['date'], df_views['count'], marker='o', linewidth=2, markersize=4, label='Total Views')
        ax1.plot(df_views['date'], df_views['uniques'], marker='s', linewidth=2, markersize=4, label='Unique Visitors')
        ax1.set_ylabel('Count', fontsize=12)
        ax1.set_title(f'Repository Views - {REPOSITORY}', fontsize=14, fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        ax1.tick_params(axis='x', rotation=45)
    
    # Clones plot
    if clones_data:
        df_clones = pd.DataFrame(clones_data)
        df_clones['date'] = pd.to_datetime(df_clones['timestamp'])
        
        ax2.plot(df_clones['date'], df_clones['count'], marker='o', linewidth=2, markersize=4, label='Total Clones', color='green')
        ax2.plot(df_clones['date'], df_clones['uniques'], marker='s', linewidth=2, markersize=4, label='Unique Cloners', color='darkgreen')
        ax2.set_ylabel('Count', fontsize=12)
        ax2.set_title('Repository Clones', fontsize=14, fontweight='bold')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        ax2.tick_params(axis='x', rotation=45)
    
    # Stars plot
    if stars_data:
        df_stars = pd.DataFrame(stars_data)
        df_stars['date'] = pd.to_datetime(df_stars['timestamp'])
        
        ax3.plot(df_stars['date'], df_stars['total'], marker='o', linewidth=2, markersize=4, label='Total Stars', color='gold')
        ax3.fill_between(df_stars['date'], df_stars['total'], alpha=0.3, color='gold')
        ax3.set_xlabel('Date', fontsize=12)
        ax3.set_ylabel('Count', fontsize=12)
        ax3.set_title('Repository Stars (Cumulative)', fontsize=14, fontweight='bold')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        ax3.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        ax3.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(TRAFFIC_DIR / 'traffic_plot.png', dpi=150, bbox_inches='tight')
    print("Plot saved to .traffic/traffic_plot.png")

def create_readme(summary, views_data, clones_data, stars_data):
    """Create a dynamic README with current stats"""
    last_updated = datetime.fromisoformat(summary['last_updated']).strftime('%B %d, %Y at %H:%M UTC')
    
    # Calculate stats
    total_views = sum(v.get('count', 0) for v in views_data)
    total_unique_views = sum(v.get('uniques', 0) for v in views_data)
    total_clones = sum(c.get('count', 0) for c in clones_data)
    total_unique_clones = sum(c.get('uniques', 0) for c in clones_data)
    current_stars = stars_data[-1]['total'] if stars_data else 0
    
    # Get latest 14-day stats
    latest_views = views_data[-14:] if len(views_data) >= 14 else views_data
    latest_clones = clones_data[-14:] if len(clones_data) >= 14 else clones_data
    
    recent_views = sum(v.get('count', 0) for v in latest_views)
    recent_unique_views = sum(v.get('uniques', 0) for v in latest_views)
    recent_clones = sum(c.get('count', 0) for c in latest_clones)
    recent_unique_clones = sum(c.get('uniques', 0) for c in latest_clones)
    
    readme_content = f"""# Repository Traffic Statistics

> Last updated: **{last_updated}**

## 📊 Traffic Overview

![Traffic Stats](traffic_plot.png)

## 📈 Current Stats

### Latest 14 Days
| Metric | Total | Unique |
|--------|-------|--------|
| 👁️ **Views** | {recent_views:,} | {recent_unique_views:,} |
| 📦 **Clones** | {recent_clones:,} | {recent_unique_clones:,} |
| ⭐ **Stars** | {current_stars:,} | - |

### All Time (Since Tracking Started)
| Metric | Total | Unique | Data Points |
|--------|-------|--------|-------------|
| 👁️ **Views** | {total_views:,} | {total_unique_views:,} | {len(views_data)} days |
| 📦 **Clones** | {total_clones:,} | {total_unique_clones:,} | {len(clones_data)} days |
| ⭐ **Stars** | {current_stars:,} | - | {len(stars_data)} days |

## 📁 Available Files

- **`traffic_plot.png`** - Visual charts (views, clones, stars)
- **`traffic_chart.txt`** - ASCII charts for terminal viewing
- **`views.json`** - Historical views data
- **`clones.json`** - Historical clones data
- **`stars.json`** - Historical stars data
- **`summary.json`** - Statistics summary

## 🔄 How This Works

This data is automatically collected every 14 days by a GitHub Action. Since GitHub only retains traffic data for 14 days, this action:
1. Fetches the current data from GitHub's API
2. Merges it with historical data
3. Builds a continuous timeline beyond the 14-day limit

Each run adds to the existing dataset, creating a complete historical record.
"""
    
    with open(TRAFFIC_DIR / 'README.md', 'w') as f:
        f.write(readme_content)
    
    print("Dynamic README saved to .traffic/README.md")


def main():
    print(f"Fetching traffic data for {REPOSITORY}...")
    
    # Fetch current traffic data
    views_response = fetch_traffic_data('views')
    clones_response = fetch_traffic_data('clones')
    
    print("Fetching stargazers data...")
    stargazers = fetch_stargazers()
    stars_data = aggregate_stars_by_date(stargazers)
    
    # Load historical data (accumulates over time)
    historical_views = load_historical_data('views.json')
    historical_clones = load_historical_data('clones.json')
    historical_stars = load_historical_data('stars.json')
    
    # Merge with new data (GitHub API returns last 14 days, we merge to build complete history)
    updated_views = merge_traffic_data(historical_views, views_response.get('views', []))
    updated_clones = merge_traffic_data(historical_clones, clones_response.get('clones', []))
    updated_stars = merge_traffic_data(historical_stars, stars_data)
    
    # Save updated data
    save_data('views.json', updated_views)
    save_data('clones.json', updated_clones)
    save_data('stars.json', updated_stars)
    
    # Create ASCII charts
    ascii_content = "# Repository Traffic Stats\n"
    ascii_content += f"Last updated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}\n"
    ascii_content += f"Repository: {REPOSITORY}\n"
    ascii_content += "\n" + create_ascii_chart(updated_views, 'views')
    ascii_content += "\n" + create_ascii_chart(updated_clones, 'clones')
    ascii_content += "\n" + create_ascii_chart(updated_stars, 'new stars')
    
    with open(TRAFFIC_DIR / 'traffic_chart.txt', 'w') as f:
        f.write(ascii_content)
    
    print("ASCII chart saved to .traffic/traffic_chart.txt")
    
    # Create plot
    create_plot(updated_views, updated_clones, updated_stars)
    
    # Create summary
    summary = {
        'last_updated': datetime.now().isoformat(),
        'total_data_points': {
            'views': len(updated_views),
            'clones': len(updated_clones),
            'stars': len(updated_stars)
        },
        'latest_stats': {
            'views': updated_views[-1] if updated_views else None,
            'clones': updated_clones[-1] if updated_clones else None,
            'stars': updated_stars[-1] if updated_stars else None
        }
    }
    
    save_data('summary.json', summary)
    
    # Create dynamic README
    create_readme(summary, updated_views, updated_clones, updated_stars)
    
    print(f"Summary saved. Views: {len(updated_views)}, Clones: {len(updated_clones)}, Stars: {len(updated_stars)}")

if __name__ == '__main__':
    main()